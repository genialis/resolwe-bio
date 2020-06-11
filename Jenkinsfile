throttle(["resolwe_bio"]) {

    node {
        // NOTE: To avoid exceeding the maximum allowed shebang lenght when calling pip due very
        // long paths of Jenkins' workspaces, we need to set a shorter Tox's working directory
        // path.
        // More info: http://tox.readthedocs.io/en/latest/example/jenkins.html#avoiding-the-path-too-long-error-with-long-shebang-lines
        def tox_workdir = "${env.HOME}/.tox-${env.BUILD_TAG}"
        // Extra arguments passed to Tox.
        def tox_extra_args = ""
        if (env.BRANCH_NAME && env.BRANCH_NAME == "master" ||
            env.CHANGE_TARGET && env.CHANGE_TARGET == "master") {
            // NOTE: If we are building the "master" branch or a pull request against the "master"
            // branch, we allow installing pre-releases with the pip command.
            tox_extra_args += "--pre"
        }

        // NOTE: Tests could hang unexpectedly and never release the Jenkins executor. Thus we set
        // a general timeout for tests' execution.
        timeout(time: 90, unit: "MINUTES") {

            try {
                stage("Checkout") {
                    // Clean up the workspace directory to ensure we will do a clean git check out.
                    deleteDir()

                    // Check out the same revision as this script is loaded from.
                    checkout scm

                    // Check if the pull request is up to date.
                    if (env.CHANGE_TARGET) {
                        git_change_target_merge_base = sh (
                            script: "git merge-base HEAD origin/${env.CHANGE_TARGET}",
                            returnStdout: true
                        ).trim()

                        git_change_target_sha = sh (
                            script: "git rev-parse origin/${env.CHANGE_TARGET}",
                            returnStdout: true
                        ).trim()

                        if (git_change_target_merge_base != git_change_target_sha) {
                            error(
                                """
                                Pull request is not up-to-date!

                                Please, rebase your pull request on top of '${env.CHANGE_TARGET}'
                                (commit: ${git_change_target_sha}).
                                """.stripIndent()
                            )
                        }
                    }
                }

                stage("Prepare") {
                    // Fetch git LFS files from remote and checkout required working tree files.
                    sh "git lfs pull"
                }

                stage("Test") {
                    // NOTE: Tox environments that don't run the ordinary tests are run first so if
                    // any of them fails, developer will get the feedback right away (rather than
                    // having to wait for all the ordinary tests to run).
                    withEnv([
                        "RESOLWE_POSTGRESQL_USER=postgres",
                        // Set database name and Redis prefix to a unique value so multiple test
                        // instances can run at the same time.
                        "RESOLWE_POSTGRESQL_NAME=${env.BUILD_TAG}",
                        "RESOLWE_MANAGER_REDIS_PREFIX=resolwe-bio.manager.${env.BUILD_TAG}",
                        "RESOLWE_DOCKER_COMMAND=sudo docker",
                        // Limit the number of parallel Django test processes.
                        "DJANGO_TEST_PROCESSES=6",
                        "TOX_WORKDIR=${tox_workdir}"
                    ]) {
                        withEnv([
                            // NOTE: These ports are set to telnet's port (23) to ensure the
                            // following Tox environments don't require access to these services.
                            "RESOLWE_POSTGRESQL_PORT=23",
                            "RESOLWE_ES_PORT=23",
                            "RESOLWE_REDIS_PORT=23"
                        ]) {
                            sh "tox -e docs ${tox_extra_args}"

                            sh "tox -e linters ${tox_extra_args}"

                            sh "tox -e packaging ${tox_extra_args}"

                            sh "tox -e extra ${tox_extra_args}"
                        }
                        withEnv([
                            // NOTE: These ports must correspond to project's services running on
                            // the Jenkins server.
                            "RESOLWE_POSTGRESQL_PORT=55433",
                            "RESOLWE_ES_PORT=59201",
                            "RESOLWE_REDIS_PORT=56380"
                        ]) {
                            sh "tox -e migrations ${tox_extra_args}"

                            if (env.CHANGE_TARGET) {
                                // Run partial test suite depending on detected changes to the
                                // target branch if we are testing a pull request.
                                withEnv([
                                    "RESOLWE_TEST_ONLY_CHANGES_TO=origin/${env.CHANGE_TARGET}"
                                ]) {
                                    sh "tox -e py36-partial ${tox_extra_args}"
                                }
                            } else {
                                sh "tox -e py36 ${tox_extra_args}"
                            }
                        }
                    }
                }

                stage("Clean") {
                    // Clean up the workspace directory after a successful build to free the disk
                    // space.
                    deleteDir()
                }

            } catch (e) {
                currentBuild.result = "FAILED"
                // Report failures only when testing the "master" or a "stable/*" branch.
                if (env.BRANCH_NAME == "master" || env.BRANCH_NAME.startsWith("stable/") ) {
                    notifyFailed()
                }
                throw e
            } finally {
                // Manually remove Tox's working directory since it is created outside Jenkins's
                // workspace.
                sh "rm -rf ${tox_workdir}"
            }
        }
    }
}

def notifyFailed() {
    slackSend(
        color: "#FF0000",
        message: "FAILED: Job ${env.JOB_NAME} (build #${env.BUILD_NUMBER}) ${env.BUILD_URL}"
    )
}
