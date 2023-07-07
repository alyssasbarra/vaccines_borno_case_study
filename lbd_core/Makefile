# targets for Jenkinsfile
# use bash and not sh (default) - https://stackoverflow.com/a/589300
SHELL := FILEPATH

# Get latest test_deploy image
newest_image := $(shell ls -t FILEPATH | head -n1)

# Define the CLUSTER variable (mocking prod-el7)
CLUSTER:= "prod-el7"

## Core Jenkinsfile/CI targets

# Unifying recipe
ci-all:	ci-verify-env ci-run-tests

# Verify environment
.PHONY: ci-verify-env
ci-verify-env:
	SINGULARITYENV_SGE_ENV=$(CLUSTER) FILEPATH -s $(newest_image) -e s "FILEPATH"

# Run R tests
.PHONY: ci-run-tests
ci-run-tests:
	SINGULARITYENV_SGE_ENV=$(CLUSTER) FILEPATH -s $(newest_image) -e s "FILEPATH"

# Create documentation
.PHONY: ci-make-docs
ci-make-docs:
	if [[ "$(GIT_BRANCH)" =~ merge ]]; then \
        	echo 'Making docs'; \
		SINGULARITYENV_SGE_ENV=$(CLUSTER) FILEPATH -s $(newest_image) -e s "FILEPATH"; \
	else \
		echo 'Not a merge; skipping ci-make-docs'; \
    	fi

## HELPER targets
# USERNAME's one-off targets
.PHONY: test
test:
	testing/run-tests
