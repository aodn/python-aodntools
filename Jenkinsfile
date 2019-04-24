pipeline {
    agent { dockerfile true }

    stages {
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