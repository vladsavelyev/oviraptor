from os.path import dirname, abspath


def package_path():
    return dirname(abspath(__file__))

