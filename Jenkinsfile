#!groovy

pipeline {
    agent none

    stages {
        stage('clean') {
            agent { label 'master' }
            steps {
                sh 'git clean -fdx'
            }
        }

        stage('container') {
            agent {
                dockerfile {
                    args '-v ${HOME}/.eggs:/home/builder/.eggs'
                    additionalBuildArgs '--build-arg BUILDER_UID=${JENKINS_UID:-9999}'
                }
            }
            environment {
                HOME = '/home/builder'
            }
            stages {
                stage('version') {
                    steps {
                        sh 'bumpversion patch'
                    }
                }
                stage('release') {
                    when { branch 'master' }
                    steps {
                        sh 'bumpversion --tag --commit release'
                    }
                }
                stage('test') {
                    steps {
                        sh 'python setup.py test'
                    }
                }
                stage('package') {
                    steps {
                        sh 'python setup.py bdist_wheel'
                    }
                }
            }
            post {
                success {
                    dir('dist/') {
                        archiveArtifacts artifacts: '*.whl', fingerprint: true, onlyIfSuccessful: true
                    }
                }
            }
        }
    }
}
