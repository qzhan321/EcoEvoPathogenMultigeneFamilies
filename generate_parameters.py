#!/usr/bin/env python
'''
    generate_parameters.py
    
    Generates parameters.hpp from parameters.hpp.template and input file with parameter
    values (either .py or .json format).
    Run `./generate_parameters.py -h' to see documentation of command-line options.
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

script_dir = os.path.dirname(__file__)

try:
    basestring = basestring
except NameError:
    basestring = str

def main():
    '''Runs everything.'''
    args = parse_arguments()
    
    try:
        os.makedirs(os.path.join(args.d, 'generated'))
    except:
        pass
    
    params = load_parameters(args.p)
    
    src_filename = os.path.join(script_dir, 'src', 'parameters.hpp.template')
    dst_filename = os.path.join(args.d, 'generated', 'parameters.hpp')
    
    print('Generating\n  {}\nfrom\n  {}\n'.format(dst_filename, src_filename))
    process_template(
        params,
        src_filename,
        dst_filename
    )

def parse_arguments():
    '''Parses command-line arguments.'''
    
    parser = argparse.ArgumentParser(
        description = '''
            Generates parameters.hpp from parameters.hpp.template and the provided
            parameters file, in either .json or .py format.
        ''',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-d', metavar = '<destination>',
        default = os.path.join(script_dir, 'build'),
        help = '''
            Destination directory for built model.
            Parameters file will go in <destination>/generated/parameters.hpp .
        '''
    )
    parser.add_argument(
        '-p', metavar = '<params-file>',
        default = os.path.join(script_dir, 'parameters-example.py'),
        help = '''
            Filename for model parameters, either a Python module (.py) or JSON (.json).
            These parameters are used to generate `parameters.hpp' from `parameters.hpp.template',
            replacing entries of the form {<variable_name>} with <variable_name> = <value>,
            where <value> is a valid C++ assignment expression for the variable.
            Integer numbers are converted into integer literals.
            Floating-point numbers are converted into floating-point literals,.
            Strings are converted into string literals.
            Lists are converted into C++ uniform-initialization lists: {<value1>, <value2>, ...},
            appropriate for types such as std::vector<double>.
            Nested lists such as std::vector<std::vector<double>> are also allowed.
            Dictionaries are not currently supported.
            Use parameter names as variable names in a Python module, or as dictionary keys in a JSON file.
        '''
    )
    return parser.parse_args()

def load_parameters(params_filename):
    '''Loads parameters from Python module or JSON file.'''
    
    base, ext = os.path.splitext(params_filename)
    if ext == '.py':
        return load_parameters_python(params_filename)
    elif ext == '.json':
        return load_parameters_json(params_filename)
    
    print('{}: unknown file type'.format(params_filename))
    sys.exit(1)

def load_parameters_python(params_filename):
    '''Loads parameters from Python module using importlib.'''
    
    old_sys_path = sys.path
    sys.path.insert(0, os.path.dirname(params_filename))
    module_name = os.path.splitext(os.path.basename(params_filename))[0]
    params = importlib.import_module(module_name)
    sys.path = old_sys_path
    
    return {name: getattr(params, name) for name in params.__dict__.keys() if not name.startswith('_')}
    
def load_parameters_json(params_filename):
    '''Loads parameters from JSON dictionary.'''
    with open(params_filename) as f:
        return json.load(f)

def process_template(params, src_filename, dst_filename):
    '''Replaces `{varname}' with `varname = value', where value is taken from params.
    
    Value is formatted appropriately using the format_value function.
    '''
    
    # Construct map: varname -> `varname = value` to use as a formatting dictionary
    param_value_map = {}
    for name, value in params.items():
        param_value_map[name] = '{} = {}'.format(name, format_value(name, value))
    
    with open(src_filename) as sf:
        with open(dst_filename, 'w') as df:
            df.write(sf.read().format(**param_value_map))

def format_value(varname, value):
    '''Formats a parameter value for substitution into a template.
    
    The Python type is used to format the value for C++.
    '''
    if isinstance(value, list):
        values = [format_value(varname, val) for val in value]
        values_formatted = '{{{}}}'.format(', '.join(values))
        return values_formatted
    elif value is False:
        return 'false'
    elif value is True:
        return 'true'
    elif isinstance(value, int):
        return str(value)
    elif isinstance(value, float):
        return json.dumps(value) # It's possible this will break in weird situations
    elif isinstance(value, basestring):
        return '{}'.format(value)
    print('Error: invalid type {} for value {} for parameter {}'.format(type(value), value, varname))
    sys.exit(1)

if __name__ == '__main__':
    main()
