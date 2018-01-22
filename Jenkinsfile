throttle(["resolwe_bio"]) {

    node {
        // NOTE: To avoid exceeding the maximum allowed shebang lenght when calling pip due very
        // long paths of Jenkins' workspaces, we need to set a shorter Tox's working directory path
        // More info: http://tox.readthedocs.io/en/latest/example/jenkins.html#avoiding-the-path-too-long-error-with-long-shebang-lines
        def tox_workdir = "${env.HOME}/.tox-${env.BUILD_TAG}"
        // extra arguments passed to Tox
        def tox_extra_args = ""
        if (env.BRANCH_NAME && env.BRANCH_NAME == "master" ||
            env.CHANGE_TARGET && env.CHANGE_TARGET == "master") {
            // NOTE: If we are building the "master" branch or a pull request against the "master"
            // branch, we allow installing pre-releases with the pip command.
            tox_extra_args += "--pre"
        }

        try {
            stage("Checkout") {
                // check out the same revision as this script is loaded from
                checkout scm
                // fetch git LFS files from remote and checkout required working tree files
                sh "git lfs pull"
            }

            stage("Test") {
                withEnv(["RESOLWE_POSTGRESQL_USER=postgres",
                         "RESOLWE_POSTGRESQL_PORT=55433",
                         "RESOLWE_ES_PORT=59201",
                         // set database name to a unique value
                         "RESOLWE_POSTGRESQL_NAME=${env.BUILD_TAG}",
                         "RESOLWE_REDIS_PORT=56380",
                         "RESOLWE_DOCKER_COMMAND=sudo docker",
                         "TEST_SUITE=resolwe_bio.tests.processes.test_enrichment.EnrichmentProcessorTestCase.test_go_enrichment",
                         // limit the number of parallel Django test processes
                         "DJANGO_TEST_PROCESSES=6",
                         "TOX_WORKDIR=${tox_workdir}"]) {
                    // run non-test Tox environments first so that if any of them fails, developer
                    // will get the feedback right away (rather than having to wait for all
                    // ordinary tests to run)
                    if (env.CHANGE_TARGET) {
                        // run partial test suite depending on detected changes to the target
                        // branch if we are testing a pull request
                        withEnv(["RESOLWE_TEST_ONLY_CHANGES_TO=origin/${env.CHANGE_TARGET}"]) {
                            sh "tox -e py34-partial ${tox_extra_args}"
                        }
                    } else {
                        sh "tox -e py34 ${tox_extra_args}"
                    }
                }
            }

        } catch (e) {
            currentBuild.result = "FAILED"
            // report failures only when testing the "master" branch
            if (env.BRANCH_NAME == "master") {
                notifyFailed()
            }
            throw e
        } finally {
            // manually remove Tox's working directory since it is created outside Jenkins's
            // workspace
            sh "rm -rf ${tox_workdir}"
        }
    }
}

def notifyFailed() {
    slackSend(
        color: "#FF0000",
        message: "FAILED: Job ${env.JOB_NAME} (build #${env.BUILD_NUMBER}) ${env.BUILD_URL}"
    )
}
