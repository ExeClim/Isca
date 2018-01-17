import copy
from jinja2 import Template

_TEMPLATE = Template("""
{%- macro fortrantrue(t) -%}
{%- if t -%}
.true.
{%- else -%}
.false.
{%- endif -%}
{%- endmacro -%}
"FMS Model results"
{% if calendar -%}
0001 1 1 0 0 0
{%- else -%}
0 0 0 0 0 0
{%- endif %}
# = output files =
# file_name, output_freq, output_units, format, time_units, long_name
{% for file in outputfiles %}
"{{ file.name }}", {{ file.freq }}, "{{ file.units }}", 1, "{{ file.time_units }}", "time",
{% endfor %}
# = diagnostic field entries =
# module_name, field_name, output_name, file_name, time_sampling, time_avg, other_opts, precision

{% for file in outputfiles %}
{% for field in file.fields -%}
"{{ field.module}}", "{{ field.name }}", "{{ field.name }}", "{{ file.name }}", "all", {{ fortrantrue(field.time_avg) }}, "none", 2,
{% endfor %}
{% endfor %}
""")



class DiagTable(object):
    def __init__(self):
        super(DiagTable, self).__init__()
        self.files = {}
        self.calendar = None

    def add_file(self, name, freq, units="hours", time_units=None):
        if time_units is None:
            time_units = units
        self.files[name] = {
            'name': name,
            'freq': freq,
            'units': units,
            'time_units': time_units,
            'fields': []
        }

    def add_field(self, module, name, time_avg=False, files=None):
        if files is None:
            files = self.files.keys()

        for fname in files:
            self.files[fname]['fields'].append({
                'module': module,
                'name': name,
                'time_avg': time_avg
                })

    def copy(self):
        d = DiagTable()
        d.files = copy.deepcopy(self.files)
        return d

    def has_calendar(self):
        if self.calendar is None or self.calendar.lower() == 'no_calendar':
            return False
        else:
            return True

    def write(self, filename):
        vars = {'calendar': self.has_calendar(), 'outputfiles': self.files.values()}
        _TEMPLATE.stream(**vars).dump(filename)

    def is_valid(self):
        return len(self.files) > 0

class DiagTableFile(object):
    def __init__(self, filename):
        self.infile = filename
        self.calendar = 'Undefined'

    def __repr__(self):
        return "DiagFile('%s')" % self.infile

    def write(self, filename):
        with open(filename, 'w') as outfile:
            with open(self.infile, 'r') as infile:
                outfile.write(infile.read())

    def is_valid(self):
        return True
