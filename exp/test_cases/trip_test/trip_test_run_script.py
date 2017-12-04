from trip_test_functions import run_all_tests

if __name__=="__main__":
    #Base commit is the earlier commit you want to compare against
    base_commit = '155661f8c7945049cbac0dcf2019bb17fe7a6a8d'
    #later commit is the newer commit you're wanting to test
    later_commit = 'ec29bf389cf5ac53b50b23c363040479a6392e52'


    #List of test cases to check
    exps_to_check = ['axisymmetric', 'bucket_model', 'frierson', 'giant_planet', 'held_suarez', 'MiMA', 'realistic_continents_fixed_sst', 'realistic_continents_variable_qflux', 'top_down_test', 'variable_co2_grey', 'variable_co2_rrtm']
        
    run_all_tests(base_commit, later_commit, exps_to_check, num_cores_to_use=4)

