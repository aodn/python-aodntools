#!groovy

pipeline {
    agent none

    stages {
        stage('container') {
            agent {
                dockerfile {
                    args '-v ${HOME}/bin:${HOME}/bin'
                    additionalBuildArgs '--build-arg BUILDER_UID=$(id -u)'
                }
            }
            stages {
                stage('clean') {
                    steps {
                        sh 'git reset --hard'
                        sh 'git clean -xffd'
                    }
                }
                stage('release') {
                    when { branch 'master' }
                    steps {
                        withCredentials([usernamePassword(credentialsId: env.GIT_CREDENTIALS_ID, passwordVariable: 'GIT_PASSWORD', usernameVariable: 'GIT_USERNAME')]) {
                            sh './bumpversion.sh'
                        }
                    }
                }
                stage('test') {
                    steps {
                        sh 'pip install --user -r test_requirements.txt'
                        sh 'pytest'
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
