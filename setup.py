from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="eden-numerical-forensics",
    version="0.1.0",
    author="André Philipsson",
    author_email="andre@eusoueden.com",
    description="Cryptographic key auditing toolkit for RSA-like moduli",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sdneycouto-jpg/eden-numerical-forensics",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Security :: Cryptography",
    ],
    python_requires=">=3.10",
    install_requires=[
        "gmpy2>=2.1.0",
    ],
    extras_require={
        "dev": ["pytest>=7.0"],
    },
)
