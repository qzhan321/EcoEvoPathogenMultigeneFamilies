#include "varmodel.hpp"

#include "StrainManager.hpp"
#include "GeneManager.hpp"
#include "PopulationManager.hpp"
#include "HostManager.hpp"
#include "InfectionManager.hpp"
#include "ImmuneHistoryManager.hpp"
#include "LocusImmunityManager.hpp"
#include "AlleleRefManager.hpp"
#include "random.hpp"
#include "util.hpp"
#include "parameters.hpp"
#include "EventQueue.hpp"

#include <unistd.h>
#include <sstream>
#include <algorithm>
#include <bitset>
#include <chrono>
using namespace std::chrono;

namespace varmodel {

#pragma mark \
*** SIMULATION STATE ***
    
    double INF = std::numeric_limits<double>::infinity();

std::mt19937_64 rng(RANDOM_SEED);

double now = 0.0;
uint64_t recomb_counter = 0;
uint64_t functional_count = 0;
uint64_t mut_count = 0;
double next_sampling_time = HOST_SAMPLING_START_YEAR * T_YEAR + HOST_SAMPLING_PERIOD[0];
uint64_t current_host_sampling_period_counter = 1;
//always output summary every 30 days;
double next_summary_time = T_BURNIN;
double next_verification_time = VERIFICATION_ON ? 0.0 : INF;
double next_checkpoint_time = SAVE_TO_CHECKPOINT ? 0.0 : INF;
//add a new type of event that introduce global mutations to the population
double next_global_mutation_time = EXPECTED_EQUILIBRIUM;
// Get starting timepoint
auto start = high_resolution_clock::now();
auto stop = high_resolution_clock::now();

uint64_t n_infections_cumulative = 0;

std::array<uint64_t, N_LOCI> n_alleles = N_ALLELES_INITIAL;
std::array<uint64_t, 2> GROUP_GENE_NUMBER = {
    (uint64_t)round(GROUP_GENE_RATIOS[0]*N_GENES_PER_STRAIN),
    (uint64_t)round(GROUP_GENE_RATIOS[1]*N_GENES_PER_STRAIN)};


std::unordered_map<uint64_t,std::array<uint64_t,2>> pairIndexMap;

double current_pop_size = 0;

StrainManager strain_manager;

GeneManager gene_manager;
std::unordered_map<
    std::array<uint64_t, N_LOCI>,
    Gene *,
    HashArray<uint64_t, N_LOCI>
> alleles_genes_map;

PopulationManager population_manager;
HostManager host_manager;
InfectionManager infection_manager;
ImmuneHistoryManager immune_history_manager;
LocusImmunityManager locus_immunity_manager;

AlleleRefManager allele_ref_manager;
std::array<std::vector<AlleleRef *>, N_LOCI> allele_refs;

sqlite3 * sample_db;

sqlite3_stmt * summary_stmt;
sqlite3_stmt * summary_alleles_stmt;

sqlite3_stmt * host_stmt;
sqlite3_stmt * strain_stmt;
sqlite3_stmt * gene_stmt;
sqlite3_stmt * allele_stmt;

sqlite3_stmt * sampled_inf_stmt;
sqlite3_stmt * sampled_imm_stmt;
sqlite3_stmt * sampled_strain_stmt;
sqlite3_stmt * sampled_gene_stmt;
sqlite3_stmt * sampled_allele_stmt;
sqlite3_stmt * extinct_gene_stmt;
sqlite3_stmt * extinct_geneAllele_stmt;
sqlite3_stmt * sampled_duration_stmt;

#pragma mark \
*** EVENT QUEUES ***
    
    double get_biting_time(Population * p) { return p->next_biting_time; }
EventQueue<Population, get_biting_time> biting_queue; 

double get_immigration_time(Population * p) { return p->next_immigration_time; }
EventQueue<Population, get_immigration_time> immigration_queue;

double get_next_immunity_loss_time(Host * h) { return h->next_immunity_loss_time; }
EventQueue<Host, get_next_immunity_loss_time> immunity_loss_queue;

double get_death_time(Host * host) { return host->death_time; }
EventQueue<Host, get_death_time> death_queue;

double get_transition_time(Infection * infection) { return infection->transition_time; }
EventQueue<Infection, get_transition_time> transition_queue;

double get_mutation_time(Infection * infection) { return infection->mutation_time; }
EventQueue<Infection, get_mutation_time> mutation_queue;

double get_recombination_time(Infection * infection) { return infection->recombination_time; }
EventQueue<Infection, get_recombination_time> recombination_queue;

double get_clearance_time(Infection * infection) { return infection->clearance_time; }
EventQueue<Infection, get_clearance_time> clearance_queue;

double get_IRS_time(Population * p) {return p->next_IRS_rate_change_time;}
EventQueue<Population, get_IRS_time> IRS_queue;

double get_MDA_time(Population * p) {return p->next_MDA_time;}
EventQueue<Population, get_MDA_time> MDA_queue;

#pragma mark \
*** Helper function declarations ***
    
    void verify_immunity_consistency();

void initialize(bool override_seed, uint64_t random_seed);
void clean_up();
void validate_and_load_parameters();

void write_summary();
void sample_hosts();
void write_host(Host * host);
void write_strain(Strain * strain, sqlite3_stmt * s_stmt, sqlite3_stmt * g_stmt, sqlite3_stmt * a_stmt);
void write_gene(Gene * gene, sqlite3_stmt * g_stmt, sqlite3_stmt * a_stmt);
void write_sampled_infection(Host * host, Infection * infection);
void write_sampled_immunity(Host * host);
void write_duration(Host * host, Infection * infection);
void load_checkpoint(bool should_load_rng_state);
void save_checkpoint();

void save_global_state_to_checkpoint(sqlite3 * db);
void load_global_state_from_checkpoint(sqlite3 * db, bool should_load_rng_state);
void load_allele_refs();
void initialize_event_queues_from_state();

std::string get_rng_as_string();
void set_rng_from_string(std::string const & rng_str);

void initialize_sample_db();
void finalize_sample_db();

void initialize_gene_pool();
double generate_functionality(GeneSource source, double similarity);
void initialize_populations();
void initialize_population(uint64_t order);
void initialize_population_events(Population * pop);
void initialize_population_hosts(Population * pop);
void initialize_population_infections(Population * pop);

Host * create_host(Population * pop, bool newborn);
void destroy_host(Host * host);

Gene * get_gene_with_alleles(std::array<uint64_t, N_LOCI> const & alleles);
Gene * get_or_create_gene(std::array<uint64_t, N_LOCI> const & alleles, GeneSource source, double is_functional, bool group, double recombRate, bool in_pool);
Strain * generate_random_strain();
double generate_recombRate(uint64_t cat);

void retain_gene(Gene * gene);
void release_gene(Gene *gene);
Strain * create_strain();
void destroy_strain(Strain * strain);
void retain_strain(Strain * strain);
void release_strain(Strain * strain);
void get_gene_pair_recomb_rates(Infection * infection);

std::array<uint64_t, 2> choose_recombine_genes(Infection * infection);
std::array<uint64_t, N_LOCI> recombine_alleles(std::array<uint64_t, N_LOCI> const & a1, std::array<uint64_t, N_LOCI> const & a2, uint64_t breakpoint);
Gene * recombine_alleles(Gene * gene1, Gene * gene2, uint64_t breakpoint, double similarity, bool create_new_allele);
bool contains_different_genes(Strain * strain);
void choose_next_gene(Infection * infection);
void trans_imm_trade_off(Gene * gene);
Strain * recombine_strains(Strain * s1, Strain * s2);
Strain * mutate_strain(Strain * strain);
Gene * mutate_gene(Gene * gene, GeneSource source, bool in_pool);
void recombine_infection(Infection * infection, bool simBiasRecomb);
double get_gene_similarity(Gene * gene1, Gene * gene2, uint64_t breakpoint);
double get_gene_similarity(Gene * gene1, Gene * gene2, double breakpoint);
std::vector<uint64_t> parents_same_alleles(Gene * gene1, Gene * gene2);
void destroy_infection(Infection * infection);

Gene * get_current_gene(Infection * infection);
Gene * draw_random_gene();
Gene * draw_random_group_gene(uint64_t group_id);


uint64_t get_immune_allele_count(Host * host);

void gain_immunity(Host * host, Gene * gene);
void update_next_immunity_loss_time(Host * host);

void lose_random_immunity(Host * host);
void lose_immunity(Host * host, AlleleRef * allele_ref);

void update_host_infection_times(Host * host);
void update_infection_times(Infection * infection);
double get_specific_immunity_level(Host * host, Gene * gene);
uint64_t get_active_infection_count(Host * host);
uint64_t get_active_ups_infection_count(Host * host, uint64_t upsGroup);
uint64_t get_liver_infection_count(Host * host);
void infect_host(Host * host, Strain * strain);
bool compareGeneFunc(Gene * gene1, Gene * gene2);
void perform_infection_transition(Infection * infection);
void clear_infection(Infection * infection);

void update_transition_time(Infection * infection, bool initial);
double draw_activation_time(Infection * infection);
double draw_deactivation_time(Infection * infection);

void update_mutation_time(Infection * infection, bool initial);
void update_recombination_time(Infection * infection, bool initial);
void update_clearance_time(Infection * infection, bool initial);

bool do_next_event(); 

void do_verification_event();
void do_sampling_event();
void do_summary_event();
void do_checkpoint_event();
void do_IRS_event();
void do_MDA_event();
uint64_t current_distinct_genes();

void do_biting_event();
Host * draw_random_source_host(Population * pop);
Host * draw_random_destination_host(Population * src_pop);
double distanceWeightFunction(double dist);
double calDistance(Population * src_pop, Population * dest_pop);
void transmit(Host * src_host, Host * dst_host);
double get_transmission_probability(Infection * infection);

void do_immigration_event();

void do_immunity_loss_event();
void do_death_event();
void do_transition_event();
void do_mutation_event();
void do_recombination_event();
void do_clearance_event(); 
void do_global_mutation_event();
void update_global_mutation_time();

void update_biting_time(Population * pop, bool initial);
void update_immigration_time(Population * pop, bool initial);
void update_MDA_time(Population * pop);
void update_biting_rate_change(Population * pop);
double draw_exponential_after_now(double lambda);
double draw_exponential(double lambda);
uint64_t draw_uniform_index(uint64_t size);
uint64_t draw_uniform_index_except(uint64_t size, uint64_t except_index);
double draw_uniform_real(double min, double max);
double draw_normal(double mean, double var, bool limit);
std::vector<uint64_t> draw_uniform_indices_without_replacement(uint64_t n, uint64_t k);
bool draw_bernoulli(double p);
uint64_t draw_discrete_distribution(std::vector<double> weights);
#pragma mark \
*** Printing/debugging helpers ***
    
    uint64_t same_time_count = 0;
uint64_t call_depth = 0;

#define BEGIN() {                              \
call_depth++;                                  \
if(PRINT_FUNCTION_TRACE) {                     \
    for(uint64_t i = 0; i < call_depth; i++) { \
        printf(" ");                           \
    }                                          \
    printf("BEGIN: %s\n", __func__);           \
}                                              \
}

#define RETURN(x) {                            \
if(PRINT_FUNCTION_TRACE) {                     \
    for(uint64_t i = 0; i < call_depth; i++) { \
        printf(" ");                           \
    }                                          \
    printf("RETURN: %s\n", __func__);          \
}                                              \
call_depth--;                                  \
return x;                                      \
}

#define PRINT_DEBUG(level, ...) {               \
if(level <= PRINT_DEBUG_LEVEL) {                \
    for(uint64_t i = 0; i <= call_depth; i++) { \
        printf(" ");                            \
    }                                           \
    printf(__VA_ARGS__);                        \
    printf("\n");                               \
}                                               \
}

#pragma mark \
*** Initialization function implementations ***
    
    void validate_and_load_parameters() {
        BEGIN();
        
        assert(RANDOM_SEED > 0);
        
        assert(T_YEAR > 0.0);
        assert(T_END >= 0.0);
        assert(T_BURNIN >= 0.0);
        assert(T_BURNIN <= T_END);
        //assert(T_BURNIN + 360 <= EXPECTED_EQUILIBRIUM);
        
        assert(SAMPLE_DB_FILENAME == "" || !file_exists(SAMPLE_DB_FILENAME));
        for (uint64_t i = 0; i < HOST_SAMPLING_PERIOD.size(); i++) {
            assert(HOST_SAMPLING_PERIOD[i] >= 0.0 && HOST_SAMPLING_PERIOD[i] <= T_YEAR);
        }
        
        assert(!SAVE_TO_CHECKPOINT || !file_exists(CHECKPOINT_SAVE_FILENAME));
        assert(!SAVE_TO_CHECKPOINT || CHECKPOINT_SAVE_PERIOD > 0.0);
        
        assert(!LOAD_FROM_CHECKPOINT || file_exists(CHECKPOINT_LOAD_FILENAME));
        
        assert(!VERIFICATION_ON || VERIFICATION_PERIOD > 0.0);
        
        assert(N_GENES_INITIAL >= 1);
        assert(N_GENES_PER_STRAIN >= 2);
        uint64_t k = 0;
        
        for(uint64_t i = 0; i < (N_GENES_PER_STRAIN-1); i++) {
            for(uint64_t j = i+1; j < N_GENES_PER_STRAIN; j++) {
                pairIndexMap[k] = {i,j};
                k+=1;
            }
        }
        
        assert(N_LOCI >= 1);
        assert(N_ALLELES_INITIAL.size() == N_LOCI);
        for(auto value : N_ALLELES_INITIAL) {
            assert(value >= 1);
        }
        
        assert(GENE_TRANSMISSIBILITY >= 0.0 && GENE_TRANSMISSIBILITY <= 1.0);
        assert(IMMUNITY_LOSS_RATE >= 0.0);
        
        assert(MUTATION_RATE >= 0.0);
        assert(T_LIVER_STAGE >= 0.0);
        
        assert(P_ECTOPIC_RECOMBINATION_IS_CONVERSION >= 0.0 && P_ECTOPIC_RECOMBINATION_IS_CONVERSION <= 1.0);
        assert(ECTOPIC_RECOMBINATION_RATE.size()<=2);
        assert(MEAN_HOST_LIFETIME > 0.0);
        assert(MAX_HOST_LIFETIME > 0.0);
        
        assert(MINIMUM_FUNCTION >= 0.0 && MINIMUM_FUNCTION <= 1.0);
        for(auto value : FUNCTION_CATEGORY){
            assert(value<=1 && value >0);
        }
        assert(TRANSITION_RATE_IMMUNE>TRANSITION_RATE_NOT_IMMUNE);
        assert(RHO>=0.01 && RHO<=0.9);
        
        assert(N_POPULATIONS >= 1);
        for(auto value : N_HOSTS) {
            assert(value >= 1);
        }
        
        assert(BITING_RATE_MEAN.size() == N_POPULATIONS);
        for(auto value : BITING_RATE_MEAN) {
            assert(value >= 0.0);
        }
        assert(BITING_RATE_RELATIVE_AMPLITUDE.size() == N_POPULATIONS);
        for(auto value : BITING_RATE_RELATIVE_AMPLITUDE) {
            assert(value >= 0.0 && value <= 1.0);
        }
        assert(BITING_RATE_PEAK_PHASE.size() == N_POPULATIONS);
        for(auto value : BITING_RATE_PEAK_PHASE) {
            assert(value >= 0.0 && value <= 1.0);
        }
        
        assert(!IMMIGRATION_ON || IMMIGRATION_RATE.size() == N_POPULATIONS);
        if(IMMIGRATION_ON) {
            for(auto value : IMMIGRATION_RATE) {
                assert(value >= 0.0);
            }
        }
        
        assert(DISTANCE_MAT.size() == N_POPULATIONS);
        for (auto distRow : DISTANCE_MAT){
            assert(distRow.size() == N_POPULATIONS);
        }
        /*if(SELECTION_MODE == GENERAL_IMMUNITY) {
         assert(GENERAL_IMMUNITY_PARAMS[0] > 0.0);
         assert(GENERAL_IMMUNITY_PARAMS[1] > 0.0);
         assert(GENERAL_IMMUNITY_PARAMS[2] > 0.0);
         assert(GENERAL_IMMUNITY_PARAMS[3] > 0.0);
        }*/
        
        RETURN();
    }

void initialize(bool override_seed, uint64_t random_seed) {
    BEGIN();
    validate_and_load_parameters();
    initialize_sample_db();
    if(override_seed) {
        rng.seed(random_seed);
    }
    
    if(LOAD_FROM_CHECKPOINT) {
        //do not use the seed from the previous run!!
        load_checkpoint(override_seed);
    }
    else {
        initialize_gene_pool();
        initialize_populations();
    }
    
    RETURN();
}

void clean_up() {
    BEGIN();
    
    finalize_sample_db();
    
    RETURN();
}

void initialize_gene_pool() {
    BEGIN();
    
    // Create AlleleRefs
    for(uint64_t i = 0; i < N_LOCI; i++) {
        assert(allele_refs[i].size() == 0);
        for(uint64_t j = 0; j < n_alleles[i]; j++) {
            AlleleRef * allele_ref = allele_ref_manager.create();
            allele_ref->locus = i;
            allele_ref->allele = j;
            allele_ref->originTime = now;
            allele_refs[i].push_back(allele_ref);
        }
    }
    
    // Create gene pool
    
    for(uint64_t i = 0; i < N_GENES_INITIAL; i++) {
        // Draw random alleles different from all other existing genes
        std::array<uint64_t, N_LOCI> alleles;
        uint64_t group = draw_discrete_distribution(GROUP_GENE_RATIOS);
        //if groups don't share alleles, they have separate sets of alleles
        if(GROUPS_DO_NOT_SHARE_ALLELE){
            do {
                for(uint64_t j = 0; j < N_LOCI; j++) {
                    alleles[j] = draw_uniform_index(n_alleles[j]/(GROUP_GENE_RATIOS[0]+GROUP_GENE_RATIOS[1])*GROUP_GENE_RATIOS[group])+group*n_alleles[j]/(GROUP_GENE_RATIOS[0]+GROUP_GENE_RATIOS[1])*GROUP_GENE_RATIOS[1-group];
                }
            } while(get_gene_with_alleles(alleles) != NULL);
            //in this implementation, we have genes that are uniform random function, and only one group of ups group
            get_or_create_gene(alleles, SOURCE_POOL_ORIGINAL, generate_functionality(SOURCE_POOL_ORIGINAL, 0), group, generate_recombRate(group),true); //some genes are less functional
            
        }else{
            do {
                for(uint64_t j = 0; j < N_LOCI; j++) {
                    alleles[j] = draw_uniform_index(n_alleles[j]);
                }
            } while(get_gene_with_alleles(alleles) != NULL);
            //in this implementation, we have genes that are uniform random function, and only one group of ups group
            get_or_create_gene(alleles, SOURCE_POOL_ORIGINAL, generate_functionality(SOURCE_POOL_ORIGINAL, 0), group, generate_recombRate(group),true); //some genes are less functional
        }
    }
    
    RETURN();
}

double generate_functionality(GeneSource source, double similarity) {
    BEGIN();
    if (source == SOURCE_RECOMBINATION) {
        if (CONTINUOUS_FUNCTION) {
            if(draw_bernoulli(similarity)){
                RETURN(draw_uniform_real(MINIMUM_FUNCTION, 1));
            }else{
                RETURN(draw_uniform_real(0, MINIMUM_FUNCTION));
            }
        }else{
            if(draw_bernoulli(similarity)){
                RETURN(1);
            }else{
                RETURN(0);
            }
        }
    }else{
        if (CONTINUOUS_FUNCTION) {
            RETURN(draw_uniform_real(0, 1));
        }else{
            RETURN(1);
        }
    }
}

void initialize_populations() {
    BEGIN();
    
    for(uint64_t i = 0; i < N_POPULATIONS; i++) {
        initialize_population(i);
    }
    
    RETURN();
}

void initialize_population(uint64_t index) {
    BEGIN();
    
    Population * pop = population_manager.create();
    pop->ind = index;
    pop->transmission_count = 0;
    pop->IRS_biting_rate = -1;
    pop->IRS_immigration_rate_factor = 1;
    pop->MDA_id = 0;
    pop->MDA_effective_period = false;
    pop->MDA_immigration_rate_factor = 1;
    pop->current_IRS_id = 0;
    pop->within_IRS_id = 0;
    pop->n_bites_cumulative = 0;
    pop->n_infected_bites = 0;
    pop->n_infected_bites_with_space = 0;
    pop->FOI = 0;
    pop->infected_ratio = 1.0;
    pop->current_pop_size = 0;
    if(IRS_ON){
        pop->next_IRS_rate_change_time = IRS_START_TIMES[pop->current_IRS_id];
        IRS_queue.add(pop);
    }
    if(MDA_ON){
        pop->next_MDA_time = MDA_START_TIMES[pop->MDA_id];
        MDA_queue.add(pop);
    }
    update_biting_time(pop, true);
    
    if(IMMIGRATION_ON) {
        update_immigration_time(pop, true);
    } else {
        pop->next_immigration_time = INF;
    }
    
    initialize_population_hosts(pop);
    initialize_population_infections(pop);
    
    RETURN();
}

void initialize_population_hosts(Population * pop) {
    BEGIN();
    
    for(uint64_t i = 0; i < N_HOSTS[pop->ind]; i++) {
        create_host(pop,false);
    }
    
    RETURN();
}

void initialize_population_infections(Population * pop) {
    BEGIN();
    
    for(uint64_t i = 0; i < N_INITIAL_INFECTIONS[pop->ind]; i++) {
        Host * host = pop->hosts.object_at_index(draw_uniform_index(pop->hosts.size())); 
        Strain * strain = generate_random_strain();
        infect_host(host, strain);
    }
    
    RETURN();
}

#pragma mark \
*** Sample database output ***
    
    void initialize_sample_db() {
        sqlite3_open(SAMPLE_DB_FILENAME.c_str(), &sample_db);
        
        sqlite3_exec(sample_db, "BEGIN TRANSACTION", NULL, NULL, NULL);
        
        sqlite3_exec(sample_db,
                     "CREATE TABLE IF NOT EXISTS hosts (id INTEGER PRIMARY KEY, population_id INTEGER, birth_time REAL, death_time REAL);",
                     NULL, NULL, NULL
        );
        
        sqlite3_exec(sample_db,
                     "CREATE TABLE IF NOT EXISTS strains (id INTEGER, ind INTEGER, gene_id INTEGER, PRIMARY KEY (id, ind));",
                     NULL, NULL, NULL
        );
        
        sqlite3_exec(sample_db,
                     "CREATE TABLE IF NOT EXISTS genes (id INTEGER PRIMARY KEY, source INTEGER, is_functional REAL, upsGroup INTEGER);",
                     NULL, NULL, NULL
        );
        
        sqlite3_exec(sample_db,
                     "CREATE TABLE IF NOT EXISTS alleles (gene_id INTEGER, locus INTEGER, allele INTEGER, PRIMARY KEY (gene_id, locus));",
                     NULL, NULL, NULL
        );
        
        sqlite3_exec(sample_db,
                     "CREATE TABLE IF NOT EXISTS summary ("
                     "time REAL, pop_id INTEGER, n_infections INTEGER, n_infected INTEGER, n_infected_bites INTEGER, n_total_bites INTEGER, n_circulating_strains INTEGER, n_circulating_genes INTEGER, exec_time REAL, n_infected_bites_with_space INTEGER, FOI INTEGER, group1NumGenes INTEGER, group2NumGenes INTEGER"
                     ");",
                     NULL, NULL, NULL
        );
        
        sqlite3_exec(sample_db,
                     "CREATE TABLE IF NOT EXISTS summary_alleles ("
                     "time REAL, pop_id INTEGER, locus INTEGER, n_circulating_alleles INTEGER"
                     ");",
                     NULL, NULL, NULL
        );
        
        sqlite3_exec(sample_db,
                     "CREATE TABLE IF NOT EXISTS sampled_infections (time REAL, host_id INTEGER, pop_id INTEGER, infection_id INTEGER, strain_id INTEGER, gene_id INTEGER, infected_length REAL, PRIMARY KEY (time, host_id, infection_id));",
                     NULL, NULL, NULL
        );
        
        sqlite3_exec(sample_db,
                     "CREATE TABLE IF NOT EXISTS sampled_immunity (time REAL, host_id INTEGER, locus INTEGER, allele INTEGER, immunity_level REAL, PRIMARY KEY (time, host_id, locus, allele));",
                     NULL, NULL, NULL
        );
        
        sqlite3_exec(sample_db,
                     "CREATE TABLE IF NOT EXISTS sampled_strains (id INTEGER, ind INTEGER, gene_id INTEGER, PRIMARY KEY (id, ind));",
                     NULL, NULL, NULL
        );
        
        sqlite3_exec(sample_db,
                     "CREATE TABLE IF NOT EXISTS sampled_genes (id INTEGER PRIMARY KEY, source INTEGER, is_functional REAL, upsGroup INTEGER, originTime REAL, recombRate REAL, deathTime REAL);",
                     NULL, NULL, NULL
        );
        
        sqlite3_exec(sample_db,
                     "CREATE TABLE IF NOT EXISTS sampled_alleles (gene_id INTEGER, locus INTEGER, allele INTEGER, originTime REAL, PRIMARY KEY (gene_id, locus));",
                     NULL, NULL, NULL
        );
        
        sqlite3_exec(sample_db,
                     "CREATE TABLE IF NOT EXISTS sampled_duration (time REAL, duration REAL, host_id INTEGER, pop_id INTEGER, infection_id INTEGER);",
                     NULL, NULL, NULL
        );
        
        sqlite3_exec(sample_db,
                     "CREATE TABLE IF NOT EXISTS extinct_genes (id INTEGER PRIMARY KEY, source INTEGER, is_functional REAL, upsGroup INTEGER, originTime REAL, recombRate REAL, deathTime REAL);",
                     NULL, NULL, NULL
        );
        
        sqlite3_exec(sample_db,
                     "CREATE TABLE IF NOT EXISTS extinct_geneAlleles (gene_id INTEGER, locus INTEGER, allele INTEGER, originTime REAL, PRIMARY KEY (gene_id, locus));",
                     NULL, NULL, NULL
        );
        
        sqlite3_prepare_v2(sample_db, "INSERT INTO summary VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?);", -1, &summary_stmt, NULL);
        sqlite3_prepare_v2(sample_db, "INSERT INTO summary_alleles VALUES (?,?,?,?);", -1, &summary_alleles_stmt, NULL);
        
        sqlite3_prepare_v2(sample_db, "INSERT INTO hosts VALUES (?,?,?,?);", -1, &host_stmt, NULL);
        sqlite3_prepare_v2(sample_db, "INSERT INTO strains VALUES (?,?,?);", -1, &strain_stmt, NULL);
        sqlite3_prepare_v2(sample_db, "INSERT INTO genes VALUES (?,?,?,?);", -1, &gene_stmt, NULL);
        sqlite3_prepare_v2(sample_db, "INSERT INTO alleles VALUES (?,?,?);", -1, &allele_stmt, NULL);
        sqlite3_prepare_v2(sample_db, "INSERT INTO sampled_infections VALUES (?,?,?,?,?,?,?);", -1, &sampled_inf_stmt, NULL);
        sqlite3_prepare_v2(sample_db, "INSERT INTO sampled_immunity VALUES (?,?,?,?,?);", -1, &sampled_imm_stmt, NULL);
        sqlite3_prepare_v2(sample_db, "INSERT OR IGNORE INTO sampled_strains VALUES (?,?,?);", -1, &sampled_strain_stmt, NULL);
        sqlite3_prepare_v2(sample_db, "INSERT OR IGNORE INTO sampled_genes VALUES (?,?,?,?,?,?,?);", -1, &sampled_gene_stmt, NULL);
        sqlite3_prepare_v2(sample_db, "INSERT OR IGNORE INTO sampled_alleles VALUES (?,?,?,?);", -1, &sampled_allele_stmt, NULL);
        sqlite3_prepare_v2(sample_db, "INSERT OR IGNORE INTO extinct_genes VALUES (?,?,?,?,?,?,?);", -1, &extinct_gene_stmt, NULL);
        sqlite3_prepare_v2(sample_db, "INSERT OR IGNORE INTO extinct_geneAlleles VALUES (?,?,?,?);", -1, &extinct_geneAllele_stmt, NULL);
        sqlite3_prepare_v2(sample_db, "INSERT OR IGNORE INTO sampled_duration VALUES (?,?,?,?,?);",-1, &sampled_duration_stmt, NULL);
        sqlite3_exec(sample_db, "COMMIT", NULL, NULL, NULL);
        sqlite3_exec(sample_db, "BEGIN TRANSACTION", NULL, NULL, NULL);
    }

void finalize_sample_db() {
    sqlite3_exec(sample_db, "COMMIT", NULL, NULL, NULL);
    
    sqlite3_finalize(summary_stmt);
    sqlite3_finalize(summary_alleles_stmt);
    
    sqlite3_finalize(host_stmt);
    sqlite3_finalize(strain_stmt);
    sqlite3_finalize(gene_stmt);
    sqlite3_finalize(sampled_gene_stmt);
    sqlite3_finalize(sampled_strain_stmt);
    sqlite3_finalize(allele_stmt);
    sqlite3_finalize(sampled_allele_stmt);
    sqlite3_finalize(extinct_gene_stmt);
    sqlite3_finalize(extinct_geneAllele_stmt);
    sqlite3_finalize(sampled_inf_stmt);
    sqlite3_finalize(sampled_imm_stmt);
    sqlite3_finalize(sampled_duration_stmt);
    sqlite3_close(sample_db);
}

#pragma mark \
*** Object management functions ***
    
    Host * create_host(Population * pop, bool newborn) {
        BEGIN();
        
        Host * host = host_manager.create();
        host->population = pop;
        host->MDA_effective_period = false;
        double lifetime = std::min(
            draw_exponential(1.0 / MEAN_HOST_LIFETIME),
            MAX_HOST_LIFETIME
        );
        if (newborn){
            host->birth_time = now;
        }else{
            host->birth_time = -draw_uniform_real(0, lifetime);
        }
        host->death_time = host->birth_time + lifetime;
        death_queue.add(host);
        
        host->completed_infection_count = 0;
        host->next_immunity_loss_time = INF;
        immunity_loss_queue.add(host);
        
        if(SELECTION_MODE == SPECIFIC_IMMUNITY) {
            ImmuneHistory * immune_history = immune_history_manager.create();
            for(uint64_t i = 0; i < N_LOCI; i++) {
                immune_history->immunity_by_locus[i] = locus_immunity_manager.create();
            }
            host->immune_history = immune_history;
        }
        else {
            host->immune_history = NULL;
        }
        
        pop->hosts.add(host);
        
        if(OUTPUT_HOSTS) {
            write_host(host);
        }
        
        RETURN(host);
    }

void destroy_host(Host * host) {
    BEGIN();
    
    Population * pop = host->population;
    PRINT_DEBUG(5, "Removing host id %llu from population %llu", host->id, pop->id);
    
    for(Infection * infection : host->infections) {
        destroy_infection(infection);
    }
    
    if(SELECTION_MODE == SPECIFIC_IMMUNITY) {
        for(LocusImmunity * immunity : host->immune_history->immunity_by_locus) {
            locus_immunity_manager.destroy(immunity);
        }
        immune_history_manager.destroy(host->immune_history);
    }
    
    pop->hosts.remove(host);
    death_queue.remove(host);
    immunity_loss_queue.remove(host);
    host_manager.destroy(host);
    
    RETURN();
}

void destroy_infection(Infection * infection) {
    BEGIN();
    release_strain(infection->strain);
    transition_queue.remove(infection);
    mutation_queue.remove(infection);
    recombination_queue.remove(infection);
    clearance_queue.remove(infection);
    infection_manager.destroy(infection);
    RETURN();
}

Gene * get_gene_with_alleles(std::array<uint64_t, N_LOCI> const & alleles) {
    BEGIN();
    
    auto itr = alleles_genes_map.find(alleles);
    if(itr == alleles_genes_map.end()) {
        RETURN(nullptr);
    }
    RETURN(itr->second);
}

Gene * get_or_create_gene(std::array<uint64_t, N_LOCI> const & alleles, GeneSource source, double is_functional, bool group, double recombRate, bool in_pool) {
    BEGIN();
    Gene * gene = get_gene_with_alleles(alleles);
    if(gene == nullptr) {
        gene = gene_manager.create();
        gene->alleles = alleles;
        gene->source = source;
        gene->is_functional = is_functional;
        gene->upsGroup = group;
        gene->originTime = now;
        gene->deathTime = -1;
        if (FUNC_DUR_DIRECT_SPECIFIFICATION) {
            gene->dur = (gene->is_functional*DURATION_CATEGORY[gene->upsGroup]);
            gene->avTransRate = (gene->is_functional*FUNCTION_CATEGORY[gene->upsGroup]);
        } else {
            trans_imm_trade_off(gene);
        }
        gene->recombRate = recombRate;
        gene->refcount = 0;
        gene->in_pool = in_pool;
        alleles_genes_map[alleles] = gene;
        if(OUTPUT_GENES) {
            write_gene(gene, gene_stmt, allele_stmt);
        }
    }
    // printf("gene new function: %f \n", gene->is_functional);
    RETURN(gene);
}

//add counting gene steps
void retain_gene(Gene * gene) {
    BEGIN();
    gene->refcount++;
    RETURN();
}

void release_gene(Gene * gene) {
    BEGIN();
    gene->refcount--;
    if (gene->refcount==0) {
        gene->deathTime = now;
        write_gene(gene, extinct_gene_stmt, extinct_geneAllele_stmt);
    }
    RETURN();
}

Strain * create_strain() {
    BEGIN();
    
    Strain * strain = strain_manager.create();
    strain->refcount = 0;
    
    RETURN(strain);
}

void destroy_strain(Strain * strain) {
    BEGIN();
    
    strain_manager.destroy(strain);
    
    RETURN();
}

void retain_strain(Strain * strain) {
    BEGIN();
    
    strain->refcount++;
    for (uint64_t i = 0; i< N_GENES_PER_STRAIN; i++){
        retain_gene(strain->genes[i]);
    }
    RETURN();
}

void release_strain(Strain * strain) {
    BEGIN();
    
    assert(strain->refcount > 0);
    strain->refcount--;
    for (uint64_t i = 0; i< N_GENES_PER_STRAIN; i++){
        release_gene(strain->genes[i]);
    }
    
    if(strain->refcount == 0)  {
        destroy_strain(strain);
    }
    
    RETURN();
}

double generate_recombRate(uint64_t cat){
    BEGIN();
    if (ECTOPIC_RECOMBINATION_RATE.size()==1) {
        RETURN(ECTOPIC_RECOMBINATION_RATE[0]);
    }else {
        if(RECOMB_UNIFORM_DIST==0){
            RETURN(draw_uniform_real(ECTOPIC_RECOMBINATION_RATE[0], ECTOPIC_RECOMBINATION_RATE[1]));
        }else if (RECOMB_UNIFORM_DIST==1){
            RETURN(ECTOPIC_RECOMBINATION_RATE[draw_uniform_index(ECTOPIC_RECOMBINATION_RATE.size())]);
        }else{
            RETURN(ECTOPIC_RECOMBINATION_RATE[cat]);
        }
    }
}

Strain * generate_random_strain() {
    BEGIN();
    
    Strain * strain = create_strain();
    auto & genes =  strain->genes;
    
    //    std::array<Gene *, N_GENES_PER_STRAIN> genes;
    
    // Old genes
    if (FIX_GENE_RATIO){
        for(uint64_t i = 0; i < GROUP_GENE_NUMBER[0]; i++) {
            genes[i] = draw_random_group_gene(0);
        }
        for(uint64_t i = GROUP_GENE_NUMBER[0]; i < N_GENES_PER_STRAIN; i++) {
            genes[i] = draw_random_group_gene(1);
        }
        
        
    }else{
        for(uint64_t i = 0; i < N_GENES_PER_STRAIN; i++) {
            genes[i] = draw_random_gene();
        }
        
    }
    
    if(OUTPUT_STRAINS) {
        write_strain(strain, strain_stmt, NULL, NULL);
    }
    
    
    RETURN(strain);
}

Strain * recombine_strains(Strain * s1, Strain * s2) {
    BEGIN();
    
    Strain * strain =  create_strain();
    auto & daughter_genes = strain->genes;
    
    // Choose random subset of all genes from strains s1 and s2
    // [0, N_GENES_PER_STRAIN) mapped to s1
    // [N_GENES_PER_STRAIN, 2 * N_GENES_PER_STRAIN) mapped to s2
    
    if(FIX_GENE_RATIO){
        // keep the number of varA and varB as a fixed ratio
        
        std::vector<Gene *> group0_genes;
        std::vector<Gene *> group1_genes;
        for(Gene * gene : s1->genes){
            if(gene->upsGroup==0){
                group0_genes.push_back(gene);
            }else{
                group1_genes.push_back(gene);
            }
        }
        for(Gene * gene : s2->genes){
            if(gene->upsGroup==0){
                group0_genes.push_back(gene);
            }else{
                group1_genes.push_back(gene);
            }
        }
        
        std::bitset<2 * N_GENES_PER_STRAIN> used;
        for(uint64_t i = 0; i < GROUP_GENE_NUMBER[0]; i++) {            // Choose random unused index across
            uint64_t src_index;
            do {
                src_index = draw_uniform_index(group0_genes.size());
            } while(used[src_index]);
            used[src_index] = true;
            daughter_genes[i] = group0_genes[src_index];
        }
        
        for(uint64_t i = GROUP_GENE_NUMBER[0]; i < N_GENES_PER_STRAIN; i++) {
            // Choose random unused index across strains
            uint64_t src_index;
            do {
                src_index = draw_uniform_index(group1_genes.size());
            } while(used[src_index+GROUP_GENE_NUMBER[0]]);
            used[src_index+GROUP_GENE_NUMBER[0]] = true;
            daughter_genes[i] = group1_genes[src_index];
        }
        
    }else{
        
        std::bitset<2 * N_GENES_PER_STRAIN> used;
        for(uint64_t i = 0; i < N_GENES_PER_STRAIN; i++) {
            // Choose random unused index across strains
            uint64_t src_index;
            do {
                src_index = draw_uniform_index(2 * N_GENES_PER_STRAIN);
            } while(used[src_index]);
            used[src_index] = true;
            if(src_index < N_GENES_PER_STRAIN) {
                daughter_genes[i] = s1->genes[src_index];
            }
            else {
                daughter_genes[i] = s2->genes[src_index - N_GENES_PER_STRAIN];
            }
        }
        
    }
    
    if(OUTPUT_STRAINS) {
        write_strain(strain, strain_stmt, NULL, NULL);
    }
    
    RETURN(strain);
}

double get_gene_similarity(Gene * gene1, Gene * gene2, uint64_t breakpoint) {
    BEGIN();
    
    double p_div = 0;
    double child_div = 0;
    double rho = 0.8; //recombination tolerance;
    double avg_mutation = 5; //average number of mutations per epitope
    for(uint64_t i = 0; i < N_LOCI; i++) {
        if(gene1->alleles[i] != gene2->alleles[i]) {
            p_div += 1;
            if(i < breakpoint) {
                child_div += 1;
            }
        }
    }
    double rho_power = child_div * avg_mutation * (p_div - child_div)* avg_mutation / (p_div * avg_mutation - 1);
    double surv_prob = pow(rho, rho_power);
    
    RETURN(surv_prob);
}

double get_gene_similarity(Gene * gene1, Gene * gene2, double breakpoint) {
    BEGIN();
    
    double p_div = 0;
    double child_div = 0;
    double avg_mutation = 5; //average number of mutations per epitope
    for(uint64_t i = 0; i < N_LOCI; i++) {
        if(gene1->alleles[i] != gene2->alleles[i]) {
            p_div += 1;
            if(i < breakpoint) {
                if(breakpoint-i>1){
                    child_div += 1;
                }else{
                    child_div += breakpoint-i;
                }
                
            }
        }
    }
    double rho_power = child_div * avg_mutation * (p_div - child_div)* avg_mutation / (p_div * avg_mutation - 1);
    double surv_prob = pow(RHO, rho_power);
    
    RETURN(surv_prob);
}

std::vector<uint64_t> parents_same_alleles(Gene * gene1, Gene * gene2){
    BEGIN();
    std::vector<uint64_t> same_allele_locations;
    for(uint64_t i = 0; i < N_LOCI; i++) {
        if(gene1->alleles[i] == gene2->alleles[i]) {
            same_allele_locations.push_back(i);
        }
    }
    RETURN(same_allele_locations);
}

Strain * mutate_strain(Strain * strain) {
    BEGIN();
    uint64_t index = draw_uniform_index(N_GENES_PER_STRAIN);
    Strain * mutated_strain = create_strain();
    mutated_strain->genes = strain->genes;
    mutated_strain->genes[index] = mutate_gene(strain->genes[index], SOURCE_MUTATION, false);
    
    if(OUTPUT_STRAINS) {
        write_strain(strain, strain_stmt, NULL, NULL);
    }
    
    RETURN(mutated_strain);
}

Gene * mutate_gene(Gene * gene, GeneSource source, bool in_pool) {
    BEGIN();
    auto alleles = gene->alleles;
    uint64_t locus = draw_uniform_index(N_LOCI); 
    alleles[locus] = n_alleles[locus];
    n_alleles[locus]++;
    
    // Add new AlleleRef for this allele
    AlleleRef * allele_ref = allele_ref_manager.create();
    allele_ref->locus = locus;
    allele_ref->allele = n_alleles[locus] - 1;
    allele_ref->originTime = now;
    allele_refs[locus].push_back(allele_ref);
    
    //keep the global pool size the same
    if(in_pool){
        gene->in_pool = false;
    }
    mut_count += 1;
    RETURN(get_or_create_gene(alleles, source, generate_functionality(SOURCE_MUTATION,0), gene->upsGroup, gene->recombRate, in_pool));
}


std::array<uint64_t, 2> choose_recombine_genes(Infection * infection){
    BEGIN();
    std::vector<double> myvector(infection->pair_recomb_rates.begin(), infection->pair_recomb_rates.end());
    uint64_t recombPair = draw_discrete_distribution(myvector);
    std::array<uint64_t, 2> recombGenes = pairIndexMap[recombPair];
    RETURN(recombGenes);
}

void choose_next_gene(Infection * infection){
    BEGIN();
    Gene * currentGene = infection->expression_order[infection->expression_index];
    uint64_t ups = currentGene->upsGroup;
    //std::cout<<get_active_ups_infection_count(infection->host, ups)-1<<" "<<get_active_ups_infection_count(infection->host, 1-ups)<<std::endl;
    if ((get_active_ups_infection_count(infection->host, ups)-1>get_active_ups_infection_count(infection->host, 1-ups)) &&
        (infection->expression_index< (N_GENES_PER_STRAIN-1))){
        
        for (uint64_t i = infection->expression_index+1; i< N_GENES_PER_STRAIN; i++){
            if(infection->expression_order[i]->upsGroup == 1-ups){
                infection->expression_order[infection->expression_index] = infection->expression_order[i];
                infection->expression_order[i] = currentGene;
                break;
            }
        }
    }
    RETURN();
}

void recombine_infection(Infection * infection, bool simBiasRecomb) {
    BEGIN();
    
    if(N_GENES_PER_STRAIN == 1) {
        RETURN();
    }
    
    Gene * new_gene_1;
    Gene * new_gene_2;
    
    // Choose genes to recombine based on their inherent rates of recombination
    std::array<uint64_t, 2> chosenRecomb = choose_recombine_genes(infection);
    uint64_t exp_index_1 = chosenRecomb[0];
    uint64_t exp_index_2 = chosenRecomb[1];
    Gene * src_gene_1 = infection->expression_order[exp_index_1];
    Gene * src_gene_2 = infection->expression_order[exp_index_2];
    if(src_gene_1 == src_gene_2) {
        RETURN();
    }
    
    bool is_conversion = draw_bernoulli(P_ECTOPIC_RECOMBINATION_IS_CONVERSION);
    
    if (ECTOPIC_RECOMBINATION_CREATE_NEW_ALLELE) {
        double breakpoint;
        if(simBiasRecomb){
            std::vector<uint64_t> same_alleles = parents_same_alleles(src_gene_1,src_gene_2);
            if(same_alleles.size()>=1){
                breakpoint = same_alleles[draw_uniform_index(same_alleles.size())]+draw_uniform_real(0, 1);
            }else{
                RETURN();
            }
        }else{
            breakpoint = draw_uniform_real(0, double(N_LOCI));
        }
        //similarity takes into account the functionality of the source genes
        double similarity = get_gene_similarity(src_gene_1, src_gene_2, breakpoint);
        // printf("similarity is %f \n", similarity);
        
        bool create_new_allele = draw_bernoulli(P_ECTOPIC_RECOMBINATION_CREATE_NEW_ALLELE);
        new_gene_1 = recombine_alleles(src_gene_1, src_gene_2, breakpoint, similarity, create_new_allele);
        // printf("gene 1 function: %f \n", new_gene_1->is_functional);
        // Under conversion, the second gene remains unchanged
        if(is_conversion) {
            new_gene_2 = src_gene_2;
        }
        // Under normal combination, both genes have material swapped
        else {
            new_gene_2 = recombine_alleles(src_gene_2, src_gene_1, breakpoint, similarity, create_new_allele);
        }
        // printf("gene 2 function: %f \n", new_gene_1->is_functional);
    } else {
        uint64_t breakpoint = draw_uniform_index(N_LOCI);
        
        // If breakpoint == 0, very little to do
        if(breakpoint == 0) {
            if(is_conversion) {
                new_gene_1 = src_gene_1;
                new_gene_2 = src_gene_1;
            }
            else {
                new_gene_1 = src_gene_1;
                new_gene_2 = src_gene_2;
            }
        } else {
            double similarity = get_gene_similarity(src_gene_1, src_gene_2, breakpoint);
            auto rec_alleles_1 = recombine_alleles(src_gene_1->alleles, src_gene_2->alleles, breakpoint);
            Gene * rec_gene_1 = get_or_create_gene(rec_alleles_1, SOURCE_RECOMBINATION, generate_functionality(SOURCE_RECOMBINATION, similarity*src_gene_1->is_functional),src_gene_1->upsGroup, src_gene_1->recombRate,false);
            new_gene_1 = rec_gene_1;
            
            // Under conversion, the second gene remains unchanged
            if(is_conversion) {
                new_gene_2 = src_gene_2;
            }
            // Under normal combination, both genes have material swapped
            else {
                auto rec_alleles_2 = recombine_alleles(src_gene_2->alleles, src_gene_1->alleles, breakpoint);
                Gene * rec_gene_2 = get_or_create_gene(rec_alleles_2, SOURCE_RECOMBINATION, generate_functionality(SOURCE_RECOMBINATION, similarity*src_gene_2->is_functional),src_gene_2->upsGroup, src_gene_2->recombRate, false);
                new_gene_2 = rec_gene_2;
            }
        }
    }
    if (new_gene_1->alleles != src_gene_1->alleles) {
        recomb_counter += 1;
        if (new_gene_1->is_functional) {
            functional_count += 1;
        }
    }
    if (new_gene_2->alleles != src_gene_2->alleles) {
        recomb_counter += 1;
        if (new_gene_2->is_functional) {
            functional_count += 1;
        }
    }
    if(RECOMBINATION_LOAD==false){
        if(new_gene_1->is_functional <= MINIMUM_FUNCTION){
            new_gene_1 = src_gene_1;
        }
        if(new_gene_2->is_functional <= MINIMUM_FUNCTION){
            new_gene_2 = src_gene_2;
        }
    }
    
    // If nothing has changed, nothing to do
    if(new_gene_1 == src_gene_1 && new_gene_2 == src_gene_2) {
        RETURN();
    }
    
    
    // Update expression order and strain
    Strain * strain = create_strain();
    auto & new_genes = strain->genes;
    new_genes = infection->strain->genes;
    if(new_gene_1 != src_gene_1) {
        infection->expression_order[exp_index_1] = new_gene_1;
        replace_first(new_genes, src_gene_1, new_gene_1);
    }
    if(new_gene_2 != src_gene_2) {
        infection->expression_order[exp_index_2] = new_gene_2;
        replace_first(new_genes, src_gene_2, new_gene_2);
    }
    Strain * old_strain  = infection->strain;
    
    infection->strain = strain;
    retain_strain(infection->strain);
    release_strain(old_strain);
    
    get_gene_pair_recomb_rates(infection);
    
    if(OUTPUT_STRAINS) {
        write_strain(strain, strain_stmt, NULL, NULL);
    }
    
    RETURN();
}

std::array<uint64_t, N_LOCI> recombine_alleles(
        std::array<uint64_t, N_LOCI> const & a1, std::array<uint64_t, N_LOCI> const & a2, uint64_t breakpoint
) {
    BEGIN();
    std::array<uint64_t, N_LOCI> arc;
    for(uint64_t i = 0; i < N_LOCI; i++) {
        if(i < breakpoint) {
            arc[i] = a1[i];
        }
        else {
            arc[i] = a2[i];
        }
    }
    RETURN(arc);
}

Gene * recombine_alleles(
        Gene * gene1, Gene * gene2, uint64_t breakpoint
    , double similarity, bool create_new_allele) {
    BEGIN();
    std::array<uint64_t, N_LOCI> & a1 = gene1->alleles;
    std::array<uint64_t, N_LOCI> & a2 = gene2->alleles;
    std::array<uint64_t, N_LOCI> arc;
    for(uint64_t i = 0; i < N_LOCI; i++) {
        if(i < breakpoint) {
            arc[i] = a1[i];
        }else if(i == breakpoint){
            //create new alleles only when the parental alleles are different
            if((create_new_allele)&&(a1[i]!=a2[i])){
                arc[i] =n_alleles[i];
                n_alleles[i]++;
                
                // Add new AlleleRef for this allele
                AlleleRef * allele_ref = allele_ref_manager.create();
                allele_ref->locus = i;
                allele_ref->allele = n_alleles[i] - 1;
                allele_ref->originTime = now;
                allele_refs[i].push_back(allele_ref);
            }else{
                arc[i] = a2[i];
            }
        }
        else {
            arc[i] = a2[i];
        }
    }
    RETURN(get_or_create_gene(arc, SOURCE_RECOMBINATION, generate_functionality(SOURCE_RECOMBINATION,similarity*gene1->is_functional), gene1->upsGroup, gene1->recombRate,false));
}

bool contains_different_genes(Strain * strain) {
    BEGIN();
    for(uint64_t i = 1; i < N_GENES_PER_STRAIN; i++) {
        if(strain->genes[i] != strain->genes[0]) {
            RETURN(true);
        }
    }
    RETURN(false);
}

Gene * draw_random_gene() {
    BEGIN();
    Gene * gene = gene_manager.objects()[draw_uniform_index(gene_manager.size())];
    while(gene->in_pool == false) {
        gene = gene_manager.objects()[draw_uniform_index(gene_manager.size())];
    }
    RETURN(gene);
}

Gene * draw_random_group_gene(uint64_t group_id) {
    BEGIN();
    Gene * gene = gene_manager.objects()[draw_uniform_index(gene_manager.size())];
    while((gene->upsGroup != group_id)||(gene->in_pool == false)){
        gene = gene_manager.objects()[draw_uniform_index(gene_manager.size())];
    }
    RETURN(gene);
}

Gene * get_current_gene(Infection * infection) {
    BEGIN();
    Gene * gene = infection->expression_order[infection->expression_index];
    RETURN(gene);
}

void trans_imm_trade_off(Gene * gene){
    BEGIN();
    //assert(gene->is_functional>=MINIMUM_FUNCTION);
    //RETURN(pow((1-pow((gene->is_functional-0.2)/0.8, 1/TRADE_OFF_S)),TRADE_OFF_S)*(15-3)+3); this tradeoff function actually doesn't work in this scenario
    //total transmission y = function*duration;
    //Total parasites = ((1+r)^nc-1)/r = 4095; r is growth rate from protein, nc is number of generations
    //duration = log((500*r+1),(1+r))*2; average transmission rate is = log((500*20+1),(1+20))*2/duration
    
    if (CONTINUOUS_FUNCTION){
        double gr = gene->is_functional*20;
        gene->dur = log(500*gr+1)/log(1+gr)*2;
        gene->avTransRate = log(500*20+1)/log(1+20)*2/gene->dur;
    }else{
        gene->dur = 1/TRANSITION_RATE_NOT_IMMUNE/(gene->is_functional*FUNCTION_CATEGORY[gene->upsGroup]);
        gene->avTransRate = (gene->is_functional*FUNCTION_CATEGORY[gene->upsGroup]);
        //var A get cleared in m days, var B/C get cleared in 2*m days
    }
    //std::cout<<gene->is_functional<<" dura "<<dur<<std::endl;
    RETURN();
}

void gain_immunity(Host * host, Gene * gene) {
    BEGIN();
    ImmuneHistory * immune_history = host->immune_history;
    //double addImm = 1+1*double(gene->is_functional);
    double addImm = 1;
    for(uint64_t i = 0; i < N_LOCI; i++) {
        LocusImmunity * immunity = immune_history->immunity_by_locus[i];
        auto itr = immunity->immunity_level_by_allele.find(gene->alleles[i]);
        if(itr == immunity->immunity_level_by_allele.end()) {
            immunity->immunity_level_by_allele[gene->alleles[i]] = addImm;
            
            AlleleRef * allele_ref = allele_refs[i][gene->alleles[i]];
            immune_history->alleles_with_immunity.add(allele_ref);
        }
        else {
            immunity->immunity_level_by_allele[gene->alleles[i]]+= addImm;
        }
    }
    update_next_immunity_loss_time(host);
    RETURN();
}

void update_next_immunity_loss_time(Host * host) {
    BEGIN();
    if(SELECTION_MODE != SPECIFIC_IMMUNITY) {
        RETURN();
    }
    //uint64_t immune_allele_count = get_immune_allele_count(host);
    host->next_immunity_loss_time = draw_exponential_after_now(IMMUNITY_LOSS_RATE * host->immune_history->alleles_with_immunity.size());
    immunity_loss_queue.update(host);
    RETURN();
}

void lose_random_immunity(Host * host) {
    BEGIN();
    
    ImmuneHistory * immune_history = host->immune_history;
    auto & alleles_with_immunity = immune_history->alleles_with_immunity; 
    assert(alleles_with_immunity.size() > 0);
    
    uint64_t index = draw_uniform_index(alleles_with_immunity.size());
    AlleleRef * ar = alleles_with_immunity.object_at_index(index);
    
    lose_immunity(host, ar);
    update_next_immunity_loss_time(host);
    
    RETURN();
}

void lose_immunity(Host * host, AlleleRef * allele_ref) {
    uint64_t locus = allele_ref->locus;
    uint64_t allele = allele_ref->allele; 
    
    ImmuneHistory * immune_history = host->immune_history;
    auto & immunity_level_by_allele = immune_history->immunity_by_locus[locus]->immunity_level_by_allele;
    uint64_t level = immunity_level_by_allele[allele];
    assert(level > 0);
    if(level == 1) {
        immunity_level_by_allele.erase(allele);
        immune_history->alleles_with_immunity.remove(allele_ref);
    }
    else {
        immunity_level_by_allele[allele]--;
    }
}

uint64_t get_immune_allele_count(Host * host) {
    BEGIN();
    RETURN(host->immune_history->alleles_with_immunity.size());
}

double get_specific_immunity_level(Host * host, Gene * gene) {
    BEGIN();
    
    double immunity_count = 0;
    for(uint64_t i = 0; i < N_LOCI; i++) {
        auto & immunity_level_by_allele = host->immune_history->immunity_by_locus[i]->immunity_level_by_allele;
        auto itr = immunity_level_by_allele.find(gene->alleles[i]);
        if(itr != immunity_level_by_allele.end()) {
            if(itr->second >= 1){
                immunity_count +=  1;
            }
        }
    }
    
    RETURN(immunity_count / (double)N_LOCI);
}

uint64_t get_active_infection_count(Host * host) {
    BEGIN();
    uint64_t count = 0;
    for(Infection * infection : host->infections) {
        if(infection->expression_index >= 0) {
            count += 1;
        }
    }
    
    RETURN(count);
}

uint64_t get_active_ups_infection_count(Host * host, uint64_t upsGroup){
    BEGIN();
    uint64_t count = 0;
    for(Infection * infection : host->infections) {
        if((infection->expression_index >= 0)&&
           (infection->expression_order[infection->expression_index]->upsGroup == upsGroup)) {
            count += 1;
        }
    }
    RETURN(count);
}

uint64_t get_liver_infection_count(Host * host) {
    BEGIN();
    uint64_t count = 0;
    for(Infection * infection : host->infections) {
        if(infection->expression_index < 0) {
            count += 1;
        }
    }
    
    RETURN(count);
}

bool compareGeneFunc(Gene * gene1, Gene * gene2){
    BEGIN();
    RETURN(gene1->avTransRate>gene2->avTransRate);
}

void get_gene_pair_recomb_rates(Infection * infection){
    BEGIN()
    
    auto & genes = infection->expression_order;
    infection->totalRecombRate = 0.0;
    double temp = 0.0;
    uint64_t k = 0;
    for(uint64_t i = 0; i < (genes.size()-1); i++) {
        for(uint64_t j = i+1; j < genes.size(); j++) {
            if(GROUPS_DO_NOT_SHARE_ALLELE){
                temp =
                    genes[i]->recombRate*genes[j]->recombRate*
                    (genes[i]->upsGroup == genes[j]->upsGroup);
            }else{
                temp =
                    genes[i]->recombRate*genes[j]->recombRate;
            }
            infection->totalRecombRate += temp;
            infection->pair_recomb_rates[k] = temp;
            k +=1;
        }
    }
    assert(k==N_GENES_PER_STRAIN*(N_GENES_PER_STRAIN-1)/2);
    // infection's total recombination rates are determined by whether in this implementation new recombination is allowed.
    infection->totalRecombRate *= RECOMBINATION_ON;
    // printf("totalRecombRate is %f \n", infection->totalRecombRate);
    RETURN()
}

void infect_host(Host * host, Strain * strain) {
    BEGIN();
    
    Infection * infection = infection_manager.create();
    infection->host = host;
    infection->strain = strain;
    retain_strain(strain);
    
    //in this implementation, genes with higher functionality expresses first.
    for(uint64_t i = 0; i < strain->genes.size(); i++) {
        infection->expression_order[i] = strain->genes[i];
    }
    
    if(HIGH_FUNCTION_EXPRESS_FIRST){
        std::sort(infection->expression_order.begin(), infection->expression_order.end(), compareGeneFunc);
    }else{
        //for each infection, the order of infection is reset
        std::shuffle(infection->expression_order.begin(), infection->expression_order.end(), rng);
    }
    
    get_gene_pair_recomb_rates(infection);
    
    if(T_LIVER_STAGE == 0.0) {
        infection->expression_index = 0;
        update_transition_time(infection, true);
    }
    else {
        //printf("transition queue size before %llu \n", transition_queue.size());
        infection->expression_index = -1;
        infection->transition_time = now + T_LIVER_STAGE;
        //printf("infection transition time before %f \n", infection->transition_time);
        transition_queue.add(infection);
        //printf("transition queue size after %llu \n", transition_queue.size());
    }
    
    update_mutation_time(infection, true);
    update_recombination_time(infection, true);
    
    if(SELECTION_MODE == GENERAL_IMMUNITY) {
        update_clearance_time(infection, true);
    }else if(SELECTION_MODE == SPECIFIC_IMMUNITY) {
        infection->clearance_time = INF;
        clearance_queue.add(infection);
    }
    
    host->infections.insert(infection);
    infection->infected_time = now;
    n_infections_cumulative++;
    
    RETURN();
}

void perform_infection_transition(Infection * infection) {
    BEGIN();
    //printf("\n            infection_id %llu:\n", infection->id);
    //printf("\n            s:");
    //printf("expression_index %llu:\n", infection->expression_index);
    
    Host * host = infection->host;
    //printf("\n            host_id %llu:\n", host->id);
    if(infection->expression_index == -1) {
        //if MOI<=4, it will start expressing, otherwise
        uint64_t actCount;
        if(UPSGROUP_DIFF_RECEPTOR){
            actCount = get_active_ups_infection_count(host, infection->expression_order[0]->upsGroup);
        }else{
            actCount = get_active_infection_count(host);
        }
        //printf("\n            actCount %llu:\n", actCount);
        //if upsgroup partition in receptor space
        //infection can start express if actCount<3 or if actCount>3,
        // the probability of start expression decreases as actCount increases
        if (UPSGROUP_DIFF_RECEPTOR) {
            if((actCount<(5-UPSGROUP_DIFF_RECEPTOR*2))||draw_bernoulli(exp(-(double(actCount)-4+UPSGROUP_DIFF_RECEPTOR*2)/2))) {
                infection->expression_index = 0;
                host->population->FOI++;
                update_mutation_time(infection, false);
                update_recombination_time(infection, false);
            }else{
                host->infections.erase(infection);
                destroy_infection(infection);
                //infection->transition_time = INF;
                //infection->expression_index = -2;//-2 denotes dormant status
                //transition_queue.update(infection);
                RETURN();
            }
        } else {
            if (actCount < MAX_ACTIVE_MOI) {
                //printf("\n            xixi:");
                infection->expression_index = 0;
                host->population->FOI++;
                update_mutation_time(infection, false);
                update_recombination_time(infection, false);
                
            }else{
                //printf("\n            haha:");
                host->infections.erase(infection);
                destroy_infection(infection);
                RETURN();
            }
        }
        
    }
    else {
        assert(
            infection->expression_index >= 0 &&
                infection->expression_index < N_GENES_PER_STRAIN
        );
        Gene * gene = get_current_gene(infection);
        
        if(SELECTION_MODE == SPECIFIC_IMMUNITY) {
            gain_immunity(host, gene);
        }
        
        infection->expression_index++;
        //in this implementation, genes under 0.2 will not be able to grow properly in the system, therefore not expressed.
        if(SELECTION_MODE == SPECIFIC_IMMUNITY) {
            while(infection->expression_index < N_GENES_PER_STRAIN){
                if(infection->expression_order[infection->expression_index]->is_functional<=MINIMUM_FUNCTION){
                    infection->expression_index++;
                }else{
                    break;
                }
            }
            if((UPSGROUP_DIFF_RECEPTOR)&&(infection->expression_index < N_GENES_PER_STRAIN)){
                choose_next_gene(infection);
                //clear infection if there are more than 4 strains occupying certain receptors
                if(get_active_ups_infection_count(infection->host, infection->expression_order[infection->expression_index]->upsGroup)>4){
                    infection->expression_index = N_GENES_PER_STRAIN;
                }
            }
        }
    }
    if(infection->expression_index == N_GENES_PER_STRAIN) {
        clear_infection(infection);
    }
    else {
        update_host_infection_times(infection->host);
    }
    //printf("\n            e:");
    //printf("expression_index update %llu:\n", infection->expression_index);
    RETURN();
}

void clear_infection(Infection * infection) {
    BEGIN();
    
    Host * host = infection->host;
    //record infection duration, for every 100 infections.
    if (n_infections_cumulative%TOTAL_CLEARED_INFECTIONS_UNIT == 0) {
        write_duration(host, infection);
    }
    host->infections.erase(infection);
    
    destroy_infection(infection);
    host->completed_infection_count++;
    
    
    if(SELECTION_MODE == GENERAL_IMMUNITY) {
        update_host_infection_times(host);
    }
    
    RETURN();
}


void update_host_infection_times(Host * host) {
    BEGIN();
    for(Infection * infection : host->infections) {
        update_infection_times(infection);
    }
    RETURN();
}

void update_infection_times(Infection * infection) {
    BEGIN();
    update_transition_time(infection, false);
    RETURN();
}

void update_transition_time(Infection * infection, bool initial) {
    BEGIN();
    if(infection->expression_index < 0) {
        RETURN();
    }
    
    assert(infection->expression_index < N_GENES_PER_STRAIN);
    
    Host * host = infection->host;
    double rate;
    if(SELECTION_MODE == SPECIFIC_IMMUNITY) {
        while((get_current_gene(infection)->is_functional<=MINIMUM_FUNCTION)&&(infection->expression_index < N_GENES_PER_STRAIN)){
            infection->expression_index++;
            if(infection->expression_index == N_GENES_PER_STRAIN) {
                infection->clearance_time = now;
                clearance_queue.update(infection);
                RETURN();
            }
        }
        Gene * active_gene = get_current_gene(infection);
        double immunity_level =  get_specific_immunity_level(host, active_gene);
        assert(immunity_level >= 0.0 && immunity_level <= 1.0);
        
        double gene_specific_rate_immune = 1/(active_gene->dur*TRADE_OFF_S+(1-TRADE_OFF_S)*(1/TRANSITION_RATE_NOT_IMMUNE));
        //in this implementation, we substitute TRANSITION_RATE_NOT_IMMUNE into function related expression time.
        //from the function trans_imm_trade_off
        rate = TRANSITION_RATE_IMMUNE * gene_specific_rate_immune / (
            TRANSITION_RATE_IMMUNE * (1.0 - immunity_level) +
                gene_specific_rate_immune * immunity_level
        );
        
        //if(COINFECTION_REDUCES_TRANSMISSION){
        //rate *= pow(get_active_infection_count(infection->host),-1); //concurrent infection increases transition time
        //}
        //a small chance that the infection got cleared, proportional to the gene's immunity
        //the clearance rate is higher for genes that higher higher specific immunity
        //if (draw_bernoulli((immunity_level)*0.05)) {
        //    infection->clearance_time = draw_exponential_after_now(rate);
        //}else{
        //    infection->clearance_time = INF;
        //}
    }
    else {
        rate = TRANSITION_RATE_NOT_IMMUNE;
    }
    // printf("rate is %f \n", rate);
    infection->transition_time = draw_exponential_after_now(rate);
    
    if(initial) {
        transition_queue.add(infection);
        clearance_queue.add(infection);
    }
    else {
        transition_queue.update(infection);
        clearance_queue.update(infection);
    }
    RETURN();
}

void update_mutation_time(Infection * infection, bool initial) {
    BEGIN();
    if(infection->expression_index>-1){
        infection->mutation_time = draw_exponential_after_now(
            MUTATION_RATE * N_GENES_PER_STRAIN * N_LOCI
        );
    }else{
        infection->mutation_time = INF;
    }
    if(initial) {
        mutation_queue.add(infection);
    }
    else {
        mutation_queue.update(infection);
    }
    RETURN();
}

void update_recombination_time(Infection * infection, bool initial) {
    BEGIN();
    //ECTOPIC_RECOMBINATION_RATE describes the overall rate of fixation of a recombination
    //in a parasite genome lineage
    if(infection->expression_index>-1){
        infection->recombination_time = draw_exponential_after_now(
            infection->totalRecombRate
        );
    }else{
        infection->recombination_time = INF;
    }
    
    if(initial) {
        recombination_queue.add(infection);
    }
    else {
        recombination_queue.update(infection);
    }
    RETURN();
}

void update_clearance_time(Infection * infection, bool initial) {
    BEGIN();
    
    assert(SELECTION_MODE == GENERAL_IMMUNITY);
    
    // TODO: replace with Qixin's function based on vector of parameters
    double rate;
    Host * host = infection->host;
    double a = GENERAL_IMMUNITY_PARAMS[0];
    double b = GENERAL_IMMUNITY_PARAMS[1];
    double c = GENERAL_IMMUNITY_PARAMS[2];
    double d = GENERAL_IMMUNITY_PARAMS[3];
    
    if(host->completed_infection_count < N_INFECTIONS_FOR_GENERAL_IMMUNITY) {
        rate = 1.0 / (a +
            b * exp(
                    -c * double(host->completed_infection_count
                    )) /
                        pow(
                            d * double(host->completed_infection_count) + 1.0,
                            d
                        )
        );
    }
    else {
        rate = CLEARANCE_RATE_IMMUNE;
    }
    
    infection->clearance_time = draw_exponential_after_now(rate)+T_LIVER_STAGE;
    if(initial) {
        clearance_queue.add(infection);
    }
    else {
        clearance_queue.update(infection);
    }
    RETURN();
}


#pragma mark \
*** Main simulation loop ***
    
    enum class EventType {
        NONE,
        VERIFICATION,
        HOST_SAMPLING,
        WRITE_SUMMARY,
        CHECKPOINT,
        BITING,
        IRS,
        MDA,
        IMMIGRATION,
        IMMUNITY_LOSS,
        DEATH,
        TRANSITION,
        MUTATION,
        RECOMBINATION,
        CLEARANCE,
        GLOBAL_MUTATE
    };

void run(bool override_seed, uint64_t random_seed) {
    BEGIN();
    // Get starting timepoint
    initialize(override_seed, random_seed);
    while(do_next_event()) {  }
    clean_up();
    
    RETURN();
}

bool do_next_event() {
    BEGIN();
    double next_event_time = INF;
    EventType next_event_type = EventType::NONE;
    // printf("transition queue next time %f \n", transition_queue.next_time());
    // Find the event type with the lowest-timed event (simple linear search)
    // This is the least elegant part of the reimplementation, but it should
    // be worth it for simplicity and reduced memory usage.
    // If this section is a speed bottleneck, re-evaluate the choice.
    // It shouldn't come even close.
    if(next_sampling_time < next_event_time) {
        next_event_time = next_sampling_time;
        next_event_type = EventType::HOST_SAMPLING;
    }
    if(next_summary_time < next_event_time) {
        next_event_time = next_summary_time;
        next_event_type = EventType::WRITE_SUMMARY;
    }
    if(VERIFICATION_ON && next_verification_time < next_event_time) {
        next_event_time = next_verification_time;
        next_event_type = EventType::VERIFICATION;
    }
    if(SAVE_TO_CHECKPOINT && next_checkpoint_time < next_event_time) {
        next_event_time = next_checkpoint_time;
        next_event_type = EventType::CHECKPOINT;
    }
    if(biting_queue.next_time() < next_event_time) {
        next_event_time = biting_queue.next_time();
        next_event_type = EventType::BITING; 
    }
    if(IRS_queue.next_time() < next_event_time) {
        next_event_time = IRS_queue.next_time();
        next_event_type = EventType::IRS;
    }
    if(MDA_queue.next_time() < next_event_time) {
        next_event_time = MDA_queue.next_time();
        next_event_type = EventType::MDA;
    }
    if(immigration_queue.next_time() < next_event_time) {
        next_event_time = immigration_queue.next_time();
        next_event_type = EventType::IMMIGRATION;
    }
    if(immunity_loss_queue.next_time() < next_event_time) {
        next_event_time = immunity_loss_queue.next_time();
        next_event_type = EventType::IMMUNITY_LOSS; 
    }
    if(death_queue.next_time() < next_event_time) {
        next_event_time = death_queue.next_time();
        next_event_type = EventType::DEATH;
    }
    if(transition_queue.next_time() < next_event_time) {
        next_event_time = transition_queue.next_time();
        next_event_type = EventType::TRANSITION;
    }
    if(mutation_queue.next_time() < next_event_time) {
        next_event_time = mutation_queue.next_time();
        next_event_type = EventType::MUTATION;
    }
    if(recombination_queue.next_time() < next_event_time) {
        next_event_time = recombination_queue.next_time();
        next_event_type = EventType::RECOMBINATION;
    }
    if(
        clearance_queue.next_time() < next_event_time
    ) {
        next_event_time = clearance_queue.next_time();
        next_event_type = EventType::CLEARANCE;
    }
    if(next_global_mutation_time < next_event_time) {
        next_event_time = next_global_mutation_time;
        next_event_type = EventType::GLOBAL_MUTATE;
    }
    
    
    PRINT_DEBUG(1, "next_event_time: %f", next_event_time);
    PRINT_DEBUG(1, "next_event_type: %d", next_event_type);
    
    // Execute the next event unless it's past T_END, in which case just advance time
    // if (next_event_time < now) {
    //     printf("next_event_type: %d", next_event_type);
    // }
    // printf("next_event_type: %d \n", next_event_type);
    // printf("next_event_time: %f \n", next_event_time);
    assert(next_event_time >= now);
    if(next_event_time > T_END) {
        now = T_END;
        RETURN(false);
    }
    
    now = next_event_time;
    switch(next_event_type) {
    case EventType::NONE:
        break;
    case EventType::VERIFICATION:
        do_verification_event();
        break;
    case EventType::HOST_SAMPLING:
        do_sampling_event();
        break;
    case EventType::WRITE_SUMMARY:
        do_summary_event();
        break;
    case EventType::CHECKPOINT:
        do_checkpoint_event();
        break;
    case EventType::BITING:
        do_biting_event();
        break;
    case EventType::IRS:
        do_IRS_event();
        break;
    case EventType::MDA:
        do_MDA_event();
        break;
    case EventType::IMMIGRATION:
        do_immigration_event();
        break;
    case EventType::IMMUNITY_LOSS:
        do_immunity_loss_event();
        break;
    case EventType::DEATH:
        do_death_event();
        break;
    case EventType::TRANSITION:
        do_transition_event();
        break;
    case EventType::MUTATION:
        do_mutation_event();
        break;
    case EventType::RECOMBINATION:
        do_recombination_event();
        break;
    case EventType::CLEARANCE:
        do_clearance_event();
        break;
    case EventType::GLOBAL_MUTATE:
        do_global_mutation_event();
        break;
    }
    RETURN(true);
}

void do_verification_event() {
    BEGIN();
    verify_simulation_state();
    next_verification_time += VERIFICATION_PERIOD;
    RETURN();
}

void do_sampling_event() {
    BEGIN();
    sample_hosts();
    next_sampling_time = HOST_SAMPLING_START_YEAR * T_YEAR + (current_host_sampling_period_counter / HOST_SAMPLING_PERIOD.size()) * T_YEAR + HOST_SAMPLING_PERIOD[current_host_sampling_period_counter % HOST_SAMPLING_PERIOD.size()];
    current_host_sampling_period_counter++;
    RETURN();
}

void do_summary_event() {
    BEGIN();
    write_summary();
    next_summary_time += 30;
    RETURN();
}
void do_checkpoint_event() {
    BEGIN();
    save_checkpoint();
    next_checkpoint_time += CHECKPOINT_SAVE_PERIOD;
    RETURN();
}

void do_biting_event() {
    BEGIN();
    assert(biting_queue.size() > 0);
    Population * pop = biting_queue.head();
    PRINT_DEBUG(1, "biting event pop: %llu", pop->id);
    Host * src_host = draw_random_source_host(pop); 
    Host * dst_host = draw_random_destination_host(pop);
    pop->n_bites_cumulative++;
    
    //count the total number of infections per host
    //uint64_t srcInf = get_active_infection_count(src_host);
    transmit(src_host, dst_host);
    
    
    // Update biting event time
    update_biting_time(pop, false);
    
    RETURN();
}

void do_IRS_event() {
    BEGIN();
    assert(IRS_queue.size() > 0);
    Population * pop = IRS_queue.head();
    PRINT_DEBUG(1, "schedule next IRS event timing for pop %llu", pop->id);
    update_biting_rate_change(pop);
    update_biting_time(pop,false);
    if(IMMIGRATION_ON) {
        update_immigration_time(pop,false);
    }
    RETURN();
}

void do_MDA_event() {
    BEGIN();
    assert(MDA_queue.size() > 0);
    Population * pop = MDA_queue.head();
    PRINT_DEBUG(1, "schedule next MDA event timing for pop %llu", pop->id);
    update_MDA_time(pop);
    if(IMMIGRATION_ON) {
        update_immigration_time(pop,false);
    }
    RETURN();
}

void update_MDA_time(Population * pop) {
    BEGIN();
    if(pop->MDA_effective_period == false){
        //if the host is not failed, clear all the infections in the host
        for(Host * host : pop->hosts.as_vector()) {
            if(draw_bernoulli(1-HOST_FAIL_RATE[pop->MDA_id])) {
                std::vector<Infection *> InfectionsToRemove;
                for(Infection * infection : host->infections) {
                    if((infection->expression_index>-1) || (infection->transition_time < (now + DRUG_EFF_DURATION[pop->MDA_id]))) {
                        InfectionsToRemove.push_back(infection);
                        //clear_infection(infection);
                    }
                }
                for (uint64_t i = 0; i<InfectionsToRemove.size();i++){
                    clear_infection(InfectionsToRemove[i]);
                }
                //make sure that in the effective duration of MDA, host do not get infections.
                host->MDA_effective_period = true;
                PRINT_DEBUG(3, "removed infection number is %lu",InfectionsToRemove.size());
            }
        }
        pop->MDA_effective_period = true;
        pop->MDA_immigration_rate_factor = MDA_IMMIGRATION_RATE_FACTORS[pop->MDA_id];
        pop->next_MDA_time = (DRUG_EFF_DURATION[pop->MDA_id] - T_LIVER_STAGE)>0 ?  (now + (DRUG_EFF_DURATION[pop->MDA_id] - T_LIVER_STAGE)) : now;
        MDA_queue.update(pop);
    }else{//if the current state of MDA is on, then turn it off.
        for(Host * host : pop->hosts.as_vector()) {
            host->MDA_effective_period = false;
        }
        pop->MDA_effective_period = false;
        pop->MDA_immigration_rate_factor = 1;
        if (pop->MDA_id == (MDA_START_TIMES.size()-1)) {
            MDA_queue.remove(pop);
        }else{
            pop->MDA_id++;
            pop->next_MDA_time = MDA_START_TIMES[pop->MDA_id];
            MDA_queue.update(pop);
        }
    }
    RETURN();
}

void update_biting_rate_change(Population * pop){
    BEGIN();
    //if within_IRS_id has surpassed the vector size of one IRS rate change range, then go to the next IRS event
    //and set the biting rate to be the same as back to biting rate mean.
    if (pop->within_IRS_id == BITING_RATE_FACTORS[pop->current_IRS_id].size()) {
        if (pop->current_IRS_id == (IRS_START_TIMES.size()-1)) {
            IRS_queue.remove(pop);
        }else{
            pop->current_IRS_id++;
            pop->within_IRS_id = 0;
            pop->next_IRS_rate_change_time = IRS_START_TIMES[pop->current_IRS_id];
            IRS_queue.update(pop);
            
        }
        pop->IRS_biting_rate = -1;
        pop->IRS_immigration_rate_factor = 1;
    }else{
        pop->IRS_biting_rate = BITING_RATE_FACTORS[pop->current_IRS_id][pop->within_IRS_id];
        pop->IRS_immigration_rate_factor = IRS_IMMIGRATION_RATE_FACTORS[pop->current_IRS_id];
        pop->next_IRS_rate_change_time += RATE_FACTOR_SCALE;//update the rate daily
        IRS_queue.update(pop);
        pop->within_IRS_id++;
    }
    PRINT_DEBUG(2, "IRS biting time %f", now);
    PRINT_DEBUG(2, "next IRS biting time change is %f", pop->next_IRS_rate_change_time);
    RETURN();
    
}


Host * draw_random_source_host(Population * pop) {
    BEGIN();
    Host * host = pop->hosts.object_at_index(draw_uniform_index(pop->hosts.size()));
    RETURN(host);
}

double distanceWeightFunction(double dist){
    BEGIN();
    assert(dist > 0.0);
    RETURN(pow(dist, -DIST_POWER));
}

Host * draw_random_destination_host(Population * src_pop) {
    BEGIN();
    
    Host * host = nullptr;
    Population * dest_pop = nullptr;
    if (N_POPULATIONS>1){
        std::vector<double> weights;
        for (uint64_t i = 0; i<N_POPULATIONS; i++){
            Population * pop = population_manager.object_for_id(i);
            double dist = distanceWeightFunction(DISTANCE_MAT[src_pop->id][pop->id]);
            weights.push_back(
                dist
                * N_HOSTS[pop->id]
            );
        }
        dest_pop = population_manager.object_for_id(draw_discrete_distribution(weights));
    }else{
        dest_pop = population_manager.object_for_id(0);
    }
    host = dest_pop->hosts.object_at_index(draw_uniform_index(N_HOSTS[dest_pop->id]));
    assert(host != nullptr);
    
    RETURN(host);
}

void transmit(Host * src_host, Host * dst_host) {
    BEGIN();
    //count the total number of infections per host
    uint64_t srcInf = get_active_infection_count(src_host);
    uint64_t dstInf = get_liver_infection_count(dst_host);
    uint64_t remainSpace = MAX_LIVER_MOI-dstInf;
    if (srcInf>0) {
        src_host->population->n_infected_bites ++;
        //if the host is in MDA effective state, do not transmit
        //if the host has many strains in liver stage, then do not transmit
        if ((dst_host->MDA_effective_period)||(remainSpace <= 0)) {
            RETURN();
        }
    } else {
        RETURN();
    }
    src_host->population->n_infected_bites_with_space ++;
    std::vector<Strain *> src_strains;
    for(Infection * infection : src_host->infections) {
        if(infection->expression_index >= 0) {
            if(draw_bernoulli(get_transmission_probability(infection))) {
                src_strains.push_back(infection->strain);
            }
        }
    }
    
    // Form set of strains to transmit: some recombinants; some unmodified
    std::vector<Strain *> strains_to_transmit(src_strains.size());
    
    // Produce a set of strains of the same size as src_strains
    for(uint64_t i = 0; i < src_strains.size(); i++) {
        Strain * strain1 = src_strains[draw_uniform_index(src_strains.size())];
        Strain * strain2 = src_strains[draw_uniform_index(src_strains.size())];
        
        // If they're the same, use them unchanged. Otherwise, recombine.
        if(strain1 == strain2) {
            strains_to_transmit[i] = strain1;
        }
        else {
            strains_to_transmit[i] = recombine_strains(strain1, strain2);
        }
    }
    
    for(Strain * strain : strains_to_transmit) {
        infect_host(dst_host, strain);
    }
    RETURN();
}

double get_transmission_probability(Infection * infection) {
    BEGIN();
    if(COINFECTION_REDUCES_TRANSMISSION) {
        if(UPSGROUP_DIFF_RECEPTOR){
            RETURN(GENE_TRANSMISSIBILITY /
                get_active_ups_infection_count(infection->host, infection->expression_order[infection->expression_index]->upsGroup)*
                    infection->expression_order[infection->expression_index]->avTransRate
            );//gene function influences transmissibility
            
        }else{
            RETURN(
                GENE_TRANSMISSIBILITY /
                    get_active_infection_count(infection->host)*
                        infection->expression_order[infection->expression_index]->avTransRate
            );//gene function influences transmissibility
        }
    }
    //from growth rate simulations, it doesn't seem that multiple infections will impair transmissibility that much.
    RETURN(GENE_TRANSMISSIBILITY*infection->expression_order[infection->expression_index]->avTransRate);
}

void do_immigration_event() {
    BEGIN();
    assert(IMMIGRATION_ON);
    assert(immigration_queue.size() > 0);
    Population * pop = immigration_queue.head();
    
    PRINT_DEBUG(1, "immigration event pop: %llu", pop->id);
    
    Strain * strain = generate_random_strain();
    uint64_t index = draw_uniform_index(pop->hosts.size());
    Host * host = pop->hosts.object_at_index(index);
    
    if (get_liver_infection_count(host)<MAX_LIVER_MOI) {
        infect_host(host, strain);
    }
    
    // Update immigration event time
    update_immigration_time(pop, false);
    RETURN();
}

void do_global_mutation_event() {
    BEGIN();
    
    if(now>EXPECTED_EQUILIBRIUM) {
        mutate_gene(draw_random_gene(),SOURCE_POOL_MUTATION, true);
        // check if the current pool size is larger than expected
        printf("global mutation of new genes\n");
    }
    update_global_mutation_time();
    
    RETURN();
    
}

void do_immunity_loss_event() {
    BEGIN();
    assert(immunity_loss_queue.size() > 0);
    
    Host * host = immunity_loss_queue.head();
    
    switch(SELECTION_MODE) {
    case SPECIFIC_IMMUNITY:
        lose_random_immunity(host);
        break;
    case GENERAL_IMMUNITY:
    case NEUTRALITY:
        assert(false);
        break;
    }
    RETURN();
}

void do_death_event() {
    BEGIN();
    assert(death_queue.size() > 0);
    
    Host * host = death_queue.head();
    Population * pop = host->population;
    destroy_host(host);
    Host * new_host = create_host(pop,true);
    PRINT_DEBUG(1, "Created new host id %llu in population %llu", new_host->id, pop->id);
    RETURN();
}

void do_transition_event() {
    BEGIN();
    assert(transition_queue.size() > 0);
    
    Infection * infection = transition_queue.head();
    
    perform_infection_transition(infection);
    RETURN();
}

void do_mutation_event() {
    BEGIN();
    assert(mutation_queue.size() > 0);
    
    Infection * infection = mutation_queue.head();
    Strain * old_strain  = infection->strain;
    infection->strain = mutate_strain(old_strain);
    retain_strain(infection->strain);
    get_gene_pair_recomb_rates(infection);
    release_strain(old_strain);
    
    update_mutation_time(infection, false);
    update_transition_time(infection, false);
    update_next_immunity_loss_time(infection->host);
    
    RETURN();
}

void do_recombination_event() {
    BEGIN();
    assert(recombination_queue.size() > 0);
    
    Infection * infection = recombination_queue.head();
    //in this implementation, similarity bias recombination. only gene that share same allele can recombine
    recombine_infection(infection, SIMILARITY_BIAS_RECOMBINATION);
    
    update_recombination_time(infection, false);
    update_transition_time(infection, false);
    update_next_immunity_loss_time(infection->host);
    
    RETURN();
}

void do_clearance_event() {
    BEGIN();
    //assert(SELECTION_MODE == GENERAL_IMMUNITY);
    Infection * infection = clearance_queue.head();
    clear_infection(infection);
    RETURN();
}

#pragma mark \
*** Verification function implementations ***
    
    void verify_simulation_state() {
        BEGIN();
        
        assert(biting_queue.verify_heap());
        assert(immigration_queue.verify_heap());
        assert(immunity_loss_queue.verify_heap());
        assert(death_queue.verify_heap());
        assert(transition_queue.verify_heap());
        assert(mutation_queue.verify_heap());
        assert(recombination_queue.verify_heap());
        assert(clearance_queue.verify_heap());
        assert(IRS_queue.verify_heap());
        assert(MDA_queue.verify_heap());
        
        verify_immunity_consistency();
        
        RETURN();
    }

void verify_immunity_consistency() {
    BEGIN();
    
    for(Population * pop : population_manager.objects()) {
        for(Host * host : pop->hosts.as_vector()) {
            ImmuneHistory * immune_history = host->immune_history;
            
            // Ensure that every AlleleRef in alleles_with_immunity has an immunity count present
            for(AlleleRef * ar : immune_history->alleles_with_immunity.as_vector()) {
                assert(immune_history->immunity_by_locus[ar->locus]->immunity_level_by_allele[ar->allele] > 0);
            }
            
            // Ensure that the total number of present alleles is the same as the size of alleles_with_immunity,
            // and all present alleles have immunity levels > 0
            uint64_t allele_count = 0;
            for(uint64_t i = 0; i < N_LOCI; i++) {
                LocusImmunity * locus_immunity = immune_history->immunity_by_locus[i]; 
                allele_count += locus_immunity->immunity_level_by_allele.size();
                for(auto kv : locus_immunity->immunity_level_by_allele) {
                    assert(kv.second > 0);
                }
            }
            assert(allele_count == immune_history->alleles_with_immunity.size());
        }
    }
    
    RETURN();
}

#pragma mark \
*** Sampling implementations ***
    
    void sample_hosts() {
        BEGIN();
        
        std::vector<Host *> infected_hosts;
        for(Host * host : host_manager.objects()) {
            if(draw_bernoulli(double(HOST_SAMPLE_SIZE)/double(host_manager.size()))){
                write_sampled_immunity(host);
                for(Infection * infection : host->infections) {
                    if (infection->expression_index >= 0) {
                        write_sampled_infection(host, infection);
                    }
                }
            }
        }
        
        // Wrap up the transaction, including all the hosts, strains, genes created since the last sampling event
        sqlite3_exec(sample_db, "COMMIT", NULL, NULL, NULL);
        
        // Begin a new transaction to include all the hosts, strains, genes created before the next sampling event
        sqlite3_exec(sample_db, "BEGIN TRANSACTION", NULL, NULL, NULL);
        
        RETURN();
    }

void write_summary() {
    BEGIN();
    
    printf("\nSummary at t = %f:\n", now);
    printf("number of gene recombination %llu \n", recomb_counter);
    printf("number of functional gene recombination %llu \n", functional_count);
    printf("number of gene mutation %llu \n", mut_count);
    printf("    n_infections_cumulative: %llu\n", n_infections_cumulative);
    
    // Get duration. Substart timepoints to
    // get durarion. To cast it to proper unit
    // use duration cast method
    stop = high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dur_milli = stop - start;
    start = stop;
    
    double temp_cps = 0;
    for(Population * pop : population_manager.objects()) {
        uint64_t n_infected = 0;
        uint64_t n_infections = 0;
        
        std::unordered_set<Strain *> distinct_strains;
        std::unordered_set<Gene *> distinct_genes;
        std::array<std::unordered_set<uint64_t>, N_LOCI> distinct_alleles;
        
        uint64_t group1NumGenes = 0;
        uint64_t group2NumGenes = 0;
        for(Host * host : pop->hosts.as_vector()) {
            uint64_t actInfCount = get_active_infection_count(host);
            if(actInfCount > 0) {
                n_infected++;
            }
            n_infections += actInfCount;
            
            for(Infection * infection : host->infections) {
                if(infection->expression_index >= 0){
                    Strain * strain = infection->strain;
                    distinct_strains.insert(strain);
                    for(Gene * gene : strain->genes) {
                        distinct_genes.insert(gene);
                        for(uint64_t i = 0; i < N_LOCI; i++) {
                            distinct_alleles[i].insert(gene->alleles[i]);
                        }
                        if (gene->upsGroup == 0) {
                            group1NumGenes += 1;
                        } else {
                            group2NumGenes += 1;
                        }
                    }
                } 
            }
        }
        // printf("liver infs %llu \n", liver_infs_count);
        
        printf("\n            population %llu:\n", pop->id);
        printf("               n_infections: %llu\n", n_infections);
        printf("                 n_infected: %llu\n", n_infected);
        printf("         n_infectious_bites: %llu\n", pop->n_infected_bites);
        printf("      n_circulating_strains: %lu\n", distinct_strains.size());
        printf("        n_circulating_genes: %lu\n", distinct_genes.size());
        printf("      n_circulating_alleles: %lu,%lu\n", distinct_alleles[0].size(),distinct_alleles[1].size());
        printf("             execution time: %f\n", dur_milli.count());
        printf("             group1NumGenes: %llu\n", group1NumGenes);
        printf("             group2NumGenes: %llu\n", group2NumGenes);
        
        sqlite3_bind_double(summary_stmt, 1, now); // time
        sqlite3_bind_int64(summary_stmt, 2, pop->id); //id of the population
        sqlite3_bind_int64(summary_stmt, 3, n_infections); // n_infections
        sqlite3_bind_int64(summary_stmt, 4, n_infected); // n_infected
        sqlite3_bind_int64(summary_stmt, 5, pop->n_infected_bites); //number of infected bites
        sqlite3_bind_int64(summary_stmt, 6, pop->n_bites_cumulative); //number of total bites in the sampling period
        sqlite3_bind_int64(summary_stmt, 7, distinct_strains.size()); // n_circulating_strains
        sqlite3_bind_int64(summary_stmt, 8, distinct_genes.size()); // n_circulating_genes
        sqlite3_bind_double(summary_stmt, 9, dur_milli.count()); //execution time in milliseconds for this sampling period
        sqlite3_bind_int64(summary_stmt, 10, pop->n_infected_bites_with_space);
        sqlite3_bind_int64(summary_stmt, 11, pop->FOI);
        sqlite3_bind_int64(summary_stmt, 12, group1NumGenes);
        sqlite3_bind_int64(summary_stmt, 13, group2NumGenes);
        sqlite3_step(summary_stmt);
        sqlite3_reset(summary_stmt);
        
        for(uint64_t i = 0; i < N_LOCI; i++) {
            sqlite3_bind_double(summary_alleles_stmt, 1, now); // time
            sqlite3_bind_int64(summary_alleles_stmt, 2, pop->id); //id of the population
            sqlite3_bind_int64(summary_alleles_stmt, 3, i); // locus
            sqlite3_bind_int64(summary_alleles_stmt, 4, distinct_alleles[i].size()); // n_circulating_alleles
            sqlite3_step(summary_alleles_stmt);
            sqlite3_reset(summary_alleles_stmt);
        }
        
        if (pop->n_bites_cumulative>0) {
            pop->infected_ratio = (float)pop->n_infected_bites/(float)pop->n_bites_cumulative;
            //update immigration rate to incorporate infected
            if(IMMIGRATION_ON) {
                update_immigration_time(pop, false);
            }
        }
        pop->n_infected_bites = 0;
        pop->n_bites_cumulative = 0;
        pop->n_infected_bites_with_space = 0;
        pop->FOI = 0;
        //record year average pop size
        pop->current_pop_size = n_infections;
        temp_cps += n_infections;
    }
    temp_cps /= population_manager.size();
    current_pop_size = temp_cps;
    if (now>=EXPECTED_EQUILIBRIUM){
        update_global_mutation_time();
    }
    RETURN();
}

void write_host(Host * host) {
    sqlite3_bind_int64(host_stmt, 1, host->id);
    sqlite3_bind_int64(host_stmt, 2, host->population->id);
    sqlite3_bind_double(host_stmt, 3, host->birth_time);
    sqlite3_bind_double(host_stmt, 4, host->death_time);
    sqlite3_step(host_stmt);
    sqlite3_reset(host_stmt);
}

void write_sampled_infection(Host * host, Infection * infection) {
    Strain * strain = infection->strain;
    
    sqlite3_bind_double(sampled_inf_stmt, 1, now);
    sqlite3_bind_int64(sampled_inf_stmt, 2, host->id);
    sqlite3_bind_int64(sampled_inf_stmt, 3, host->population->id);
    sqlite3_bind_int64(sampled_inf_stmt, 4, infection->id);
    sqlite3_bind_int64(sampled_inf_stmt, 5, strain->id);
    if (infection->expression_index>-1){
        sqlite3_bind_int64(sampled_inf_stmt, 6, get_current_gene(infection)->id);
    }else{
        sqlite3_bind_int64(sampled_inf_stmt, 6, -1);
    }
    sqlite3_bind_double(sampled_inf_stmt, 7, (now-infection->infected_time));
    
    sqlite3_step(sampled_inf_stmt);
    sqlite3_reset(sampled_inf_stmt);
    
    write_strain(strain, sampled_strain_stmt, sampled_gene_stmt, sampled_allele_stmt);
}

void write_sampled_immunity(Host * host) {
    for(uint64_t i = 0; i < N_LOCI; i++) {
        auto & immunity_level_by_allele = host->immune_history->immunity_by_locus[i]->immunity_level_by_allele;
        for (auto itr = immunity_level_by_allele.begin(); itr != immunity_level_by_allele.end(); ++ itr){
            sqlite3_bind_double(sampled_imm_stmt, 1, now);
            sqlite3_bind_int64(sampled_imm_stmt, 2, host->id);
            sqlite3_bind_int64(sampled_imm_stmt, 3, i);
            sqlite3_bind_int64(sampled_imm_stmt, 4, itr->first);
            sqlite3_bind_double(sampled_imm_stmt, 5, itr->second);
            sqlite3_step(sampled_imm_stmt);
            sqlite3_reset(sampled_imm_stmt);
        }
    }
}


void write_strain(Strain * strain, sqlite3_stmt * s_stmt, sqlite3_stmt * g_stmt, sqlite3_stmt * a_stmt) {
    for(uint64_t i = 0; i < N_GENES_PER_STRAIN; i++) {
        Gene * gene = strain->genes[i];
        sqlite3_bind_int64(s_stmt, 1, strain->id);
        sqlite3_bind_int64(s_stmt, 2, i);
        sqlite3_bind_int64(s_stmt, 3, gene->id);
        sqlite3_step(s_stmt);
        sqlite3_reset(s_stmt);
        
        if(g_stmt != NULL) {
            write_gene(gene, g_stmt, a_stmt);
        }
    }
}

void write_gene(Gene * gene, sqlite3_stmt * g_stmt, sqlite3_stmt * a_stmt) {
    sqlite3_bind_int64(g_stmt, 1, gene->id);
    sqlite3_bind_int64(g_stmt, 2, gene->source);
    sqlite3_bind_double(g_stmt, 3, gene->is_functional);
    sqlite3_bind_int64(g_stmt, 4, gene->upsGroup);
    sqlite3_bind_double(g_stmt, 5, gene->originTime);
    sqlite3_bind_double(g_stmt, 6, gene->recombRate);
    sqlite3_bind_double(g_stmt, 7, gene->deathTime);
    
    sqlite3_step(g_stmt);
    sqlite3_reset(g_stmt);
    
    for(uint64_t j = 0; j < N_LOCI; j++) {
        sqlite3_bind_int64(a_stmt, 1, gene->id);
        sqlite3_bind_int64(a_stmt, 2, j);
        sqlite3_bind_int64(a_stmt, 3, gene->alleles[j]);
        sqlite3_bind_double(a_stmt, 4, allele_refs[j][gene->alleles[j]]->originTime);
        sqlite3_step(a_stmt);
        sqlite3_reset(a_stmt);
    }
    
}

void write_duration(Host * host, Infection * infection) {
    sqlite3_bind_double(sampled_duration_stmt, 1, now);
    sqlite3_bind_double(sampled_duration_stmt, 2, (now-infection->infected_time));
    sqlite3_bind_int64(sampled_duration_stmt, 3, host->id);
    sqlite3_bind_int64(sampled_duration_stmt, 4, host->population->id);
    sqlite3_bind_int64(sampled_duration_stmt, 5, host->completed_infection_count);
    sqlite3_step(sampled_duration_stmt);
    sqlite3_reset(sampled_duration_stmt);
    
}

#pragma mark \
*** Checkpoint function implementations ***
    
    void save_checkpoint() {
        BEGIN();
        std::string old_checkpoint_filename = CHECKPOINT_SAVE_FILENAME + "-old";
        if(file_exists(CHECKPOINT_SAVE_FILENAME)) {
            assert(!rename(CHECKPOINT_SAVE_FILENAME.c_str(), old_checkpoint_filename.c_str()));
        }
        
        sqlite3 * db;
        assert(!file_exists(CHECKPOINT_SAVE_FILENAME));
        sqlite3_open(CHECKPOINT_SAVE_FILENAME.c_str(), &db);
        sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, NULL);
        
        save_global_state_to_checkpoint(db);
        strain_manager.save_to_checkpoint(db);
        gene_manager.save_to_checkpoint(db);
        population_manager.save_to_checkpoint(db);
        host_manager.save_to_checkpoint(db);
        infection_manager.save_to_checkpoint(db);
        immune_history_manager.save_to_checkpoint(db);
        locus_immunity_manager.save_to_checkpoint(db);
        allele_ref_manager.save_to_checkpoint(db);
        
        sqlite3_exec(db, "COMMIT", NULL, NULL, NULL);
        sqlite3_close(db);
        
        if(file_exists(old_checkpoint_filename)) {
            assert(!unlink(old_checkpoint_filename.c_str()));
        }
        RETURN();
    }

void load_checkpoint(bool should_load_rng_state) {
    BEGIN();
    
    assert(file_exists(CHECKPOINT_LOAD_FILENAME));
    
    sqlite3 * db;
    sqlite3_open(CHECKPOINT_LOAD_FILENAME.c_str(), &db);
    
    load_global_state_from_checkpoint(db, should_load_rng_state);
    
    // Load all objects, minus references to other objects
    strain_manager.load_from_checkpoint(db);
    gene_manager.load_from_checkpoint(db);
    population_manager.load_from_checkpoint(db);
    host_manager.load_from_checkpoint(db);
    infection_manager.load_from_checkpoint(db);
    immune_history_manager.load_from_checkpoint(db);
    locus_immunity_manager.load_from_checkpoint(db);
    allele_ref_manager.load_from_checkpoint(db);
    
    // Resolve references to other objects
    strain_manager.resolve_references(db, gene_manager);
    gene_manager.resolve_references(db); // Does nothing unless referenes are added to Gene
    population_manager.resolve_references(db, host_manager);
    host_manager.resolve_references(db, population_manager, immune_history_manager, infection_manager);
    infection_manager.resolve_references(db, strain_manager, host_manager, gene_manager);
    immune_history_manager.resolve_references(db, locus_immunity_manager, allele_ref_manager);
    locus_immunity_manager.resolve_references(db); // Does nothing unless referenes are added to LocusImmunity
    allele_ref_manager.resolve_references(db);
    sqlite3_close(db);
    
    load_allele_refs();
    initialize_event_queues_from_state();
    
    RETURN();
}

void save_global_state_to_checkpoint(sqlite3 * db) {
    sqlite3_exec(db,
                 "CREATE TABLE global_state ("
                 "rng TEXT, "
                 "now REAL, "
                 "next_verification_time REAL, "
                 "next_checkpoint_time REAL, "
                 "next_global_mutation_time REAL, "
                 "n_infections_cumulative INTEGER,"
                 "current_pop_size REAL,"
                 "current_host_sampling_period_counter INTEGER"
                 ");",
                 NULL, NULL, NULL
    );
    
    sqlite3_stmt * stmt;
    sqlite3_prepare_v2(db, "INSERT INTO global_state VALUES (?,?,?,?,?,?,?);", -1, &stmt, NULL);
    std::string rng_str = get_rng_as_string();
    sqlite3_bind_text(stmt, 1, rng_str.c_str(), (int)(rng_str.size() + 1), SQLITE_STATIC);
    sqlite3_bind_double(stmt, 2, now);
    sqlite3_bind_double(stmt, 3, next_verification_time);
    sqlite3_bind_double(stmt, 4, next_checkpoint_time);
    sqlite3_bind_double(stmt, 5, next_global_mutation_time);
    sqlite3_bind_int64(stmt, 6, n_infections_cumulative);
    sqlite3_bind_double(stmt, 7, current_pop_size);
    sqlite3_bind_int64(stmt, 8, current_host_sampling_period_counter);
    sqlite3_step(stmt);
    sqlite3_finalize(stmt);
}


void load_global_state_from_checkpoint(sqlite3 * db, bool should_load_rng_state) {
    sqlite3_stmt * stmt;
    sqlite3_prepare_v2(db, "SELECT * FROM global_state LIMIT 1;", -1, &stmt, NULL);
    sqlite3_step(stmt);
    
    if(should_load_rng_state) {
        // Load rng from string representation
        std::string rng_str = (const char *)sqlite3_column_text(stmt, 0);
        PRINT_DEBUG(1, "rng read: %s", rng_str.c_str());
        set_rng_from_string(rng_str);
        assert(rng_str == get_rng_as_string());
    }
    
    now = sqlite3_column_double(stmt, 1);
    next_verification_time = sqlite3_column_double(stmt, 2);
    next_checkpoint_time = sqlite3_column_double(stmt, 3);
    next_global_mutation_time = sqlite3_column_double(stmt, 4);
    n_infections_cumulative = sqlite3_column_int(stmt, 5);
    current_pop_size = sqlite3_column_double(stmt, 6);
    current_host_sampling_period_counter = sqlite3_column_int(stmt, 7);
    next_sampling_time = HOST_SAMPLING_START_YEAR * T_YEAR + (current_host_sampling_period_counter / HOST_SAMPLING_PERIOD.size()) * T_YEAR + HOST_SAMPLING_PERIOD[current_host_sampling_period_counter % HOST_SAMPLING_PERIOD.size()];
    current_host_sampling_period_counter++;
    next_summary_time = now + 30;
    sqlite3_finalize(stmt);
}

void load_allele_refs() {
    BEGIN();
    
    // Put each allele into the vector for its locus
    for(AlleleRef * ar : allele_ref_manager.objects()) {
        allele_refs[ar->locus].push_back(ar);
    }
    
    // Sort each locus
    for(uint64_t i = 0; i < N_LOCI; i++) {
        std::sort(allele_refs[i].begin(), allele_refs[i].end(),
                  [](AlleleRef * o1, AlleleRef * o2) {
                      return o1->allele < o2->allele;
                  }
        );
    }
    
    // Verify everything's right
    // make n_alleles the size of allele_refs
    
    for(uint64_t i = 0; i < N_LOCI; i++) {
        n_alleles[i] = allele_refs[i].size();
        //assert(allele_refs[i].size() == n_alleles[i]);
        for(uint64_t j = 0; j < n_alleles[i]; j++) {
            assert(allele_refs[i][j]->locus == i);
            assert(allele_refs[i][j]->allele == j);
        }
    }
    
    RETURN();
}

void initialize_event_queues_from_state() {
    for(Population * pop : population_manager.objects()) {
        biting_queue.add(pop);
        immigration_queue.add(pop);
        if(IRS_ON){
            pop->next_IRS_rate_change_time = IRS_START_TIMES[pop->current_IRS_id];
            IRS_queue.add(pop);
        }
        if(MDA_ON){
            pop->next_MDA_time = MDA_START_TIMES[pop->MDA_id];
            MDA_queue.add(pop);
        }
        
    }
    
    for(Host * host : host_manager.objects()) {
        immunity_loss_queue.add(host);
        death_queue.add(host);
    }
    
    for(Infection * infection : infection_manager.objects()) {
        transition_queue.add(infection);
        mutation_queue.add(infection);
        recombination_queue.add(infection);
        clearance_queue.add(infection);
    }
}

std::string get_rng_as_string() {
    std::stringstream ss;
    ss << rng;
    return ss.str();
}

void set_rng_from_string(std::string const & rng_str) {
    std::stringstream ss;
    ss << rng_str;
    ss >> rng;
}


#pragma mark \
*** Random draw helper function implementations *** 
    
    double draw_exponential_after_now(double lambda) {
        BEGIN();
        double time;
        if(lambda == 0.0) {
            time = INF;
        }
        else {
            time = now + draw_exponential(lambda);
        }
        RETURN(time);
    }

double draw_exponential(double lambda) {
    BEGIN();
    RETURN(std::exponential_distribution<>(lambda)(rng));
}

uint64_t draw_uniform_index(uint64_t size) {
    BEGIN();
    RETURN(std::uniform_int_distribution<uint64_t>(0, size - 1)(rng));
}

uint64_t draw_uniform_index_except(uint64_t size, uint64_t except_index) {
    BEGIN();
    uint64_t index = draw_uniform_index(size - 1);
    if(index >= except_index) {
        index++;
    }
    RETURN(index);
}

double draw_uniform_real(double min, double max) {
    BEGIN();
    RETURN(std::uniform_real_distribution<double>(min, max)(rng));
}

double draw_normal(double mean, double var, bool limit){
    BEGIN();
    double n = std::normal_distribution<double>(mean, var)(rng);
    if (limit){
        n = n<1 ? n : 1;
        n = n>0 ? n : 0;
    }
    RETURN(n);
}

std::vector<uint64_t> draw_uniform_indices_without_replacement(uint64_t n, uint64_t k) {
    BEGIN();
    
    assert(k <= n);
    
    std::uniform_int_distribution<uint64_t> index_dist(0, n - 1);
    
    // If we're drawing less than half, draw indices to *include*
    std::vector<uint64_t> indices;
    if(k < n / 2) {
        std::unordered_set<uint64_t> include_indices;
        for(uint64_t i = 0; i < k; i++) {
            while(true) {
                uint64_t index = index_dist(rng);
                if(include_indices.find(index) == include_indices.end()) {
                    include_indices.insert(index);
                    break;
                }
            }
        }
        assert(include_indices.size() == k);
        indices.insert(indices.end(), include_indices.begin(), include_indices.end());
        std::sort(indices.begin(), indices.end());
    }
    // Otherwise draw indices to *exclude*
    else {
        std::unordered_set<uint64_t> exclude_indices;
        for(uint64_t i = 0; i < n - k; i++) {
            while(true) {
                uint64_t index = index_dist(rng);
                if(exclude_indices.find(index) == exclude_indices.end()) {
                    exclude_indices.insert(index);
                    break;
                }
            }
        }
        assert(exclude_indices.size() == n - k);
        for(uint64_t i = 0; i < n; i++) {
            if(exclude_indices.find(i) == exclude_indices.end()) {
                indices.push_back(i);
            }
        }
    }
    assert(indices.size() == k);
    RETURN(indices);
}

bool draw_bernoulli(double p) {
    BEGIN();
    RETURN(std::bernoulli_distribution(p)(rng));
}

uint64_t draw_discrete_distribution(std::vector<double> weights) {
    BEGIN();
    std::discrete_distribution<uint64_t> d(weights.begin(),weights.end());
    RETURN(d(rng));
}

void update_biting_time(Population * pop, bool initial) {
    BEGIN();
    
    double biting_rate = BITING_RATE_MEAN[pop->ind]*N_HOSTS[pop->ind];
    if(DAILY_BITING_RATE_DISTRIBUTION.size()==(int)T_YEAR){
        biting_rate *= DAILY_BITING_RATE_DISTRIBUTION[(int)(now)%(int)T_YEAR];
    }else{
        biting_rate *= (1.0 + BITING_RATE_RELATIVE_AMPLITUDE[pop->ind] *
            cos(2 * M_PI * ((now / T_YEAR) - BITING_RATE_PEAK_PHASE[pop->ind])));
    }
    if (pop->IRS_biting_rate>=0){
        biting_rate *= pop->IRS_biting_rate;   
    }
    pop->next_biting_time = draw_exponential_after_now(biting_rate);
    //printf("next_biting_time is: %f\n",pop->next_biting_time);
    if(initial) {
        biting_queue.add(pop);
    }
    else {
        biting_queue.update(pop);
    }
    pop->current_biting_rate = biting_rate;
    RETURN();
}

void update_immigration_time(Population * pop, bool initial) {
    BEGIN();
    //realized migration rates are adjusted to accomodate popsize*migration rate
    pop->next_immigration_time = draw_exponential_after_now(IMMIGRATION_RATE[pop->ind]*
        pop->current_biting_rate*pop->infected_ratio*
        pop->IRS_immigration_rate_factor*pop->MDA_immigration_rate_factor);
    if(initial) {
        immigration_queue.add(pop);
    }
    else {
        immigration_queue.update(pop);
    }
    RETURN();
}

void update_global_mutation_time() {
    BEGIN();
    double mutRate = ECTOPIC_RECOMBINATION_RATE[1] * ECTOPIC_RECOMBINATION_RATE[1] * N_GENES_PER_STRAIN * (N_GENES_PER_STRAIN - 1) / 2.0 * current_pop_size * P_ECTOPIC_RECOMBINATION_CREATE_NEW_ALLELE * P_GENE_INVASION * REGION_TO_LOCAL_POP_SIZE_RATIO;
    next_global_mutation_time = draw_exponential_after_now(mutRate);
    RETURN();
}

} // namespace varmodel
