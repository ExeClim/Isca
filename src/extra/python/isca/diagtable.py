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

def numorstr(x):
    """Try to parse a string into an int or float."""
    x = x.strip()
    if x.startswith('"'):
        return x.strip('"')
    try:
        ix = int(x)
        fx = float(x)
        return ix if ix == fx else fx
    except: pass
    if x.lower() == '.true.': return True
    if x.lower() == '.false.': return False
    return x

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

    @classmethod
    def from_file(cls, filename):
        lines = [l.strip() for l in open(filename)]
        lines = [l.split(',') for l in lines if not l.startswith("#")]
        dt = cls()
        dt.calendar = False
        #dt.files = [l[0] for l in lines if len(l)==7]
        with open(filename, 'r') as file:
            for line in file:
                lx = line.strip()
                if lx.startswith('#'):
                    continue
                if lx == '0001 1 1 0 0 0': 
                    dt.calendar = 'undefined'
                    continue
                ls = lx.split(',')
                vals = [numorstr(x) for x in ls]
                if len(ls) == 7:
                    dt.add_file(
                        name=vals[0],
                        freq=vals[1], 
                        units=vals[2], 
                        time_units=vals[4])
                elif len(ls) == 9:
                    dt.add_field(
                        module=vals[0], 
                        name=vals[1], 
                        time_avg=vals[5], 
                        files=[vals[3]])
        return dt