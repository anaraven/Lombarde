#!/usr/bin/python

from jinja2 import Environment
from jinja2.loaders import FileSystemLoader
import sys

env = Environment(loader = FileSystemLoader('.'),
    line_statement_prefix = "#+" #,
    # variable_start_string="$(",
    # variable_end_string=")"
    )

print env.get_template(sys.argv[1]).render()
