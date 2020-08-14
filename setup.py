import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Mathieu-Functions-Mikejmnez",
    version='0.0.1',
    author='Miguel Jimenez-Urias',
    author_email="mjimen17@jh.edu",
    description="Calculates Mathieu Functions of First kind for q real or purely imaginary",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/Mikejmnez/Mathieu-Functions",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
