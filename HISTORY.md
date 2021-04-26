# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v2_2_16] - 2021-04-26
### Notes
Added a contigs_500 in datadump for qc_score, easiest way to do it
### changed
- bifrost_assemblatron/datadump.py
## [v2_2_15] - 2021-02-17
### Notes
bifrostlib bump for datetime bug

## [v2_2_14] - 2021-02-15
### Notes
Fixed datadump on bbmap due to multispaces
### Changed
- bifrost_assemblatron/datadump.py
  

## [v2_2_9 - v2_2_13] - 2021-02-15
### Notes
Forgot to update HISTORY.md

- setup.py 
  - updated bifrostlib

## [v2_2_8] - 2021-02-11
### Notes
Switched install pip for prod to -e

## [v2_2_1] - 2020-12-17
### Notes
Changes to use the 2_1_0 schema, organizational updates, and updates to tests to make this work. Also updated the scheme for how the docker image is developed on to be from the root for local dev.

### Added
- docs/
  - history.rst
  - index.rst
  - readme.rst
- tests/
  - test_simple.py
- HISTORY.md
- setup.cfg

### Changed
- .dockerignore
- Dockerfile
- requirements.txt
- setup.py
- bifrost_assemblatron/
  - \_\_init\_\_.py
  - \_\_main\_\_.py
  - config.yaml
  - datadump.py
  - launcher.py
  - pipeline.smk
  - rule__summarize_depth.py
  - rule__summarize_variants.py
- .github/workflows
  - docker_build_and_push_to_dockerhub.yml
  - test_standard_workflow.yml -> run_tests.yml

### Removed
- test_1_standard_workflow.py
- docker-compose.dev.yaml
- docker-compose.yaml
- requirements.dev.txt