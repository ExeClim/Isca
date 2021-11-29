#!/usr/bin/env groovy

pipeline {
	agent any

	environment {
		GFDL_WORK = "${env.WORKSPACE}/_work"
		GFDL_DATA = "${env.WORKSPACE}/_data"
		GFDL_BASE = "${env.WORKSPACE}"
		GFDL_ENV = "emps-gv"
	}

	stages {
		stage('Setup') {
			steps {
				// checkout scm
				sh """
				# setup the python environment
				module load python/anaconda
				conda remove -y --all --name jenkins
				conda create -y -n jenkins python=2.7
				source activate jenkins
				cd $GFDL_BASE/src/extra/python
				pip install -r requirements.txt
				conda install -y scipy xarray tqdm pytest
				pip install -e .
				"""
			}
		}

		stage('Compile') {
			steps {
				sh """
				module load python/anaconda
				source activate jenkins
				cd $GFDL_BASE/test
				pytest --junitxml results.xml
				"""
				junit 'test/results.xml'
			}
		}

		// stage('Test') {
		// 	steps {
		// 		sh """
		// 		module load python/anaconda
		// 		source activate jenkins
		// 		cd $GFDL_BASE/exp/test_cases/trip_test
		// 		./trip_test_command_line -r ${env.GIT_URL} -e held_suarez f3dd4ec6b3587de4fe5fcdb4d61003a8f499931d ${env.GIT_COMMIT}
		// 		"""
		// 	}
		// }
	}
}