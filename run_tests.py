#!/usr/bin/env python
'''
    run_tests.py
    
    Runs all tests.
'''

import sys
if sys.version_info < (2,7):
    print('This script requires Python 2.7 or later.')
    sys.exit(1)

import os
import shutil
import importlib
import json
import argparse
import subprocess
import multiprocessing

script_dir = os.path.dirname(os.path.abspath(__file__))
print(script_dir)
print('here')
try:
    basestring = basestring
except NameError:
    basestring = str

def main():
    '''Runs everything.'''
    args = parse_arguments()
    
    test_filenames = [
        filename for filename in os.listdir(os.path.join(script_dir, 'tests'))
        if filename.endswith('.py') and filename != 'shared_test_params.py'
    ]
    print(test_filenames[0])
    dst_dirname = args.dest_dir
    
    args_list = [(dst_dirname, filename) for filename in test_filenames]
    print(args_list)
    if args.n_cores == 1:
        map(run_test, args_list)
    else:
        pool = multiprocessing.Pool(args.n_cores)
        pool.map(run_test, args_list)

def run_test(args):
    dst_dirname, test_filename = args
    
    test_name = os.path.splitext(test_filename)[0]
    print(test_name)
    test_dst_dirname = os.path.join(dst_dirname, test_name)
    os.makedirs(os.path.join(test_dst_dirname, 'output'))
    
    build(test_name, test_dst_dirname, os.path.join(script_dir, 'tests', test_filename))
    run(test_name, test_dst_dirname)

def parse_arguments():
    '''Parses command-line arguments.'''
    
    parser = argparse.ArgumentParser(
        description = '''Builds and runs each test.''',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-d', '-D', '--dest-dir', metavar = '<dest-dir>', default = 'build-tests')
    parser.add_argument('-n', '-N', '--n-cores', metavar = '<n-cores>', type = int, default = 1, help = 'Number of cores')
    
    return parser.parse_args()

def build(test_name, test_dst_dirname, test_filename):
    print('Building {}...'.format(test_name))
    
    with open(os.path.join(test_dst_dirname, 'output', 'build_output.txt'), 'w') as outfile:
        subprocess.Popen([
            os.path.join(script_dir, 'build.py'),
            '-p', test_filename,
            '-d', test_dst_dirname
        ], stdout=outfile, stderr=outfile).wait()
    
    print('...build complete.')
    
def run(test_name, test_dst_dirname):
    exec_filename = os.path.join(test_dst_dirname, 'bin', 'varRdiv')
    if os.path.exists(exec_filename):
        print('Running {}...'.format(test_name))
        
        output_dirname = os.path.join(test_dst_dirname, 'output')
        with open(os.path.join(output_dirname, 'run_output.txt'), 'w') as outfile:
            subprocess.Popen(
                [os.path.abspath(exec_filename)], cwd = output_dirname, stdout = outfile, stderr = outfile
            ).wait()
    
        print('...run complete.')
    else:
        print('Build failed; not running.')
    

if __name__ == '__main__':
    main()
