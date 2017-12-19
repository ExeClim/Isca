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
				checkout scm
				sh """
				# setup the python environment
				module load python/anaconda
				mkdir -p $WORKSPACE/.env
				conda config --append env_dirs $WORKSPACE/.env
				conda create -y -n jenkins python=2.7
				source activate jenkins
				cd $GFDL_BASE/src/extra/python
				pip install -r requirements.txt
				pip install -e .
				"""
			}
		}

		stage('Test') {
			steps {
				dir '${env.$GFDL_BASE}/exp/test_cases/trip_test'
				sh './trip_test_command_line -r ${env.GIT_URL} f6c2ced2a773865f7acfb615c105bc39a799bf20 ${env.GIT_COMMIT} held_suarez'
			}
		}
	}
}