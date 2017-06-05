throttle(["resolwe_bio"]) {

    node {

        try {
            stage("Checkout") {
                // check out the same revision as this script is loaded from
                checkout scm
                // fetch git LFS files from remote and checkout required working tree files
                sh "git lfs pull"
            }

            stage("Test") {
                // remove Tox's directory to force a rebuild of Tox's environments
                sh "rm -rf .tox/"
                withEnv(["RESOLWE_POSTGRESQL_USER=postgres",
                         "RESOLWE_POSTGRESQL_PORT=55433",
                         "RESOLWE_ES_PORT=59201",
                         // set database name to a unique value
                         "RESOLWE_POSTGRESQL_NAME=${env.BUILD_TAG}",
                         "RESOLWE_DOCKER_COMMAND=sudo docker",
                         // set number of parallel Django test processes to 6
                         "DJANGO_TEST_PROCESSES=6"]) {
                    // documentation, linters and packaging environments are run first so that if
                    // any of them fails, developer will get the feedback right away (rather than
                    // having to wait for all ordinary tests to run)
                    sh "tox -e docs"

                    sh "tox -e linters"

                    sh "tox -e packaging"

                    sh "echo 'Environment:' && python2.7 --version"
                    sh "tox -e py27"

                    sh "echo 'Environment:' && python3.4 --version"
                    sh "tox -e py34"
                }
            }

        } catch (e) {
            currentBuild.result = "FAILED"
            // report failures only when testing the "master" branch
            if (env.BRANCH_NAME == "master") {
                notifyFailed()
            }
            throw e
        }
    }
}

def notifyFailed() {
    slackSend(
        color: "#FF0000",
        message: "FAILED: Job ${env.JOB_NAME} (build #${env.BUILD_NUMBER}) ${env.BUILD_URL}"
    )
}
