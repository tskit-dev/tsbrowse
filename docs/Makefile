
# Need to set PYTHONPATH so that we pick up the local tsbrowse
PYPATH=${PWD}/..
TSBROWSE_VERSION:=$(shell PYTHONPATH=${PYPATH} \
   python3 -c 'import tsbrowse; print(tsbrowse.__version__.split("+")[0])')

dev:
	PYTHONPATH=${PYPATH} ./build.sh

dist:
	@echo Building distribution for tsbrowse version ${TSBROWSE_VERSION}
	sed -i s/__TSBROWSE_VERSION__/${TSBROWSE_VERSION}/g _config.yml
	PYTHONPATH=${PYPATH} ./build.sh

clean:
	rm -fR _build
