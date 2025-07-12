#!/usr/bin/env python
'''
    make.py
    
    Builds the model code.
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

script_dir = os.path.dirname(os.path.abspath(__file__))
#print(script_dir)

try:
    basestring = basestring
except NameError:
    basestring = str

def main():
    '''Runs everything.'''
    args = parse_arguments()
    
    build_sqlite(args.c_compiler, args.cflags)
    
    copy_sources(args.dest_dir)
    generate_managers(args.dest_dir)
    generate_parameters(args.dest_dir, args.params_file)
    build(args.dest_dir, args.cpp_compiler, args.cppflags)

def parse_arguments():
    '''Parses command-line arguments.'''
    
    parser = argparse.ArgumentParser(
        description = '''Builds the model in the specified directory using provided parameters.
            
            Copies source code into the `src' subdirectory, with 
        ''',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-p', '-P', '--params-file', metavar = '<params-file>',
        default = os.path.join(script_dir, 'parameters-example.py'),
        help = 'Path to parameters file in .py or .json format.'
    )
    parser.add_argument(
        '-d', '-D', '--dest-dir', metavar = '<dest-dir>',
        default = os.path.join(script_dir, 'build'),
        help = 'Path to destination directory for built model.'
    )
    parser.add_argument('-c', '--c-compiler', metavar = '<c-compiler>', default = 'cc', help = 'C compiler.')
    parser.add_argument('-C', '--cpp-compiler', metavar = '<compiler>', default = 'c++', help = 'C++ compiler.')
    parser.add_argument('-f', '--cflags', metavar = '<c-flags>', default = '-O2 -g ', help = 'C compiler flags.')
    parser.add_argument('-F', '--cppflags', metavar = '<cpp-flags>', default = '-O2 -g ', help = 'C++ compiler flags.')
    
    return parser.parse_args()

def build_sqlite(c_compiler, cflags):
    sqlite3_dir = os.path.join(script_dir, 'sqlite3')
    if os.path.exists(os.path.join(sqlite3_dir, 'sqlite3.o')):
        print('sqlite3 already built.')
    else:
        print('Building sqlite3...')
        subprocess.Popen(
            '{} {} -c -o sqlite3.o sqlite3.c'.format(c_compiler, cflags),
            cwd = sqlite3_dir, shell = True
        ).wait()
        print('sqlite3 build complete.')

def copy_sources(dst_dirname):
    try:
        os.makedirs(dst_dirname)
    except:
        pass
    
    shutil.copytree(os.path.join(script_dir, 'src'), os.path.join(dst_dirname, 'src'))
    
    git_dirname = os.path.join(dst_dirname, 'git')
    try:
        os.makedirs(git_dirname)
    except:
        pass
    with open(os.path.join(git_dirname, 'commit.txt'), 'w') as f:
        subprocess.Popen(
            ['git', 'rev-parse', 'HEAD'],
            stdout=f
        )
    with open(os.path.join(git_dirname, 'log.txt'), 'w') as f:
        subprocess.Popen(
            ['git', 'log'],
            stdout=f
        )
    with open(os.path.join(git_dirname, 'diff.txt'), 'w') as f:
        subprocess.Popen(
            ['git', 'diff'],
            stdout=f
        )

def generate_managers(dst_dirname):
    subprocess.Popen([
        os.path.join(script_dir, 'generate_managers.py'),
        '-d',
        dst_dirname
    ]).wait()

def generate_parameters(dst_dirname, params_filename):
    subprocess.Popen([
        os.path.join(script_dir, 'generate_parameters.py'),
        '-d', dst_dirname,
        '-p', params_filename
    ]).wait()

def build(dst_dirname, compiler_cmd, compiler_flags):
    os.makedirs(os.path.join(dst_dirname, 'bin'))
    
    sqlite3_dir = os.path.abspath(os.path.join(script_dir, 'sqlite3'))
    boost_dir = "/usr/local/include/"
    compile_cmd =  compiler_cmd + \
         ' ' + compiler_flags + \
         ' -std=c++11' + \
         ' -o bin/varRdiv' + \
         ' ' + '-ldl' + \
         ' ' + '-lpthread' + \
         ' ' + '-I ' + sqlite3_dir + \
         ' ' + os.path.join(sqlite3_dir, 'sqlite3.o') + \
         ' generated/managers/*.cpp' + \
         ' src/*.cpp' + \
         ' src/util/*.cpp' + \
         ' -I src' + \
         ' -I src/datamodel' + \
         ' -I src/managers' + \
         ' -I src/util' + \
         ' -I generated' + \
         ' -I generated/managers' + \
         ' -I ' + boost_dir
    print(compile_cmd)
    subprocess.Popen(compile_cmd, cwd = dst_dirname, shell = True).wait()

if __name__ == '__main__':
    main()
