# Downloading and cleaning Gaia EDR3 data

This code is based on [...]

### Pre-requisites
Requires Python 3.7+ and the libraries listed in requirements.txt, which can be installed with ```pip3 install -r requirements.txt``` or manually with pip/conda. A user and login must already be registered at [Gaia-ESA CAS](https://cas.cosmos.esa.int/cas/login?service=https%3A%2F%2Ftools.cosmos.esa.int%2Fprivacy%2Findex.php), and must be inserted in the file my_credentials.txt.

### Inputs
- ngc variable (line 34): Object name or coordinate
- Radius in degrees

### Description of the outputs
- Raw catalogue with all the columns from Gaia EDR3
- Catalogue applying all the corrections and quality flags (Gaia EDR3 papers)

<!-- <img width="1608" alt="HW77-new-example" src="https://user-images.githubusercontent.com/29663898/153299914-9427f073-aa49-4f6d-97ff-84bc53ec3cfe.png"> -->


## References (TBD)
- Lindegren+2021
- Riello+2021
- Cantat-Gaudin & Brandt (2010)
- Fabricius+2021
- Gaia Collaboration+2021
- Oliveira+2022
