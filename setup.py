
from setuptools import setup, find_packages  

setup(  
    name = 'ERgene',  
    version = '1.1.1',
    # keywords = ('chinesename',),  
    description = 'A python module that could analysis the ERGs',  
    license = 'MIT License',  
    install_requires = ['upsetplot'],  
    packages = ['ERgene'],  # 要打包的项目文件夹
    include_package_data=True,   # 自动打包文件夹内所有数据
    author = 'ZehuaZeng',  
    author_email = 'Starlitnightly@163.com',
    url = 'https://github.com/Starlitnightly/ERgene',
    # packages = find_packages(include=("*"),),  
)  
