(install)=

# Installing tsbrowse
The recommended way to install tsbrowse is using a Python virtual environment:
```
python3 -m venv env
source env/bin/activate
python3 -m pip install tsbrowse
```
Alternatively, for a containerised setup, you can use the pre-built [Docker](https://www.docker.com) image:
```
# Pull the tsbrowse image from Docker Hub
docker pull savitakartik/tsbrowse:0.0.3
# Optionally tag the image with a local name
docker tag savitakartik/tsbrowse:0.0.3 tsbrowse:0.0.3
```
For an introduction to tsbrowse, see [Introduction](intro.md)
