from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="chemistry-telegram-bot",
    version="1.0.0",
    author="Chemistry Education Bot Team",
    author_email="contact@chemistrybot.com",
    description="A comprehensive chemistry education Telegram bot with IUPAC conversion, periodic table, and reaction animations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/chemistry-telegram-bot",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Education",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Education",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "chemistry-bot=run:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["assets/*.png", "*.md", "*.txt"],
    },
)