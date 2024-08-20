import os

import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

def get_field(n):
    base = {
            "atmos_mod": "sphum_age",
            "longname": "sphum times",
            "units": "sec (kg/kg)",
            "numerical_representation": "grid",
            "hole_filling": "off",
            "advect_vert": "finite_volume_parabolic",
            "robert_filter": "on",
            "profile_type": ["fixed", "surface_value=0.0"]
        }
    base["atmos_mod"] += f"_{n+1}"
    base["longname"] += f" {n+1}-th moment "
    return base

def write_ft(name,n_moments):
    base_dir = os.path.dirname(os.path.realpath(__file__))
    # a CodeBase can be a directory on the computer,
    # useful for iterative development
    add_dir = "/src/extra/model/isca/"
    field_table_dir = GFDL_BASE + add_dir

    sphum ={
            "atmos_mod": "sphum",
            "longname": "specific humidity",
            "units": "kg/kg",
            "numerical_representation": "grid",
            "hole_filling": "off",
            "advect_vert": "finite_volume_parabolic",
            "robert_filter": "on",
            "profile_type": ["fixed", "surface_value=0.0"]
            }
        
    
    # Define the file name without an extension
    output_file = field_table_dir + name

    # Open the file in write mode
    with open(output_file, 'w') as file:
        # Write sphum tracer
        file.write(f'"TRACER",')
        for key, value in sphum.items():
            if type(value) == str:
                file.write(f'"{key}", "{value}"\n')
            elif type(value) == list:
                file.write(f'"{key}", "{value[0]}", "{value[1]}"\n')
        file.write('/ \n')

         
        for ind in range(n_moments):
            field = get_field(ind)
            file.write(f'"TRACER",')
            for key, value in field.items():
                if type(value) == str:
                    file.write(f'"{key}", "{value}"\n')
                elif type(value) == list:
                    file.write(f'"{key}", "{value[0]}", "{value[1]}"\n')
            file.write('/ \n')

        

if __name__ == "__main__":
    write_ft("field_table_test",1)