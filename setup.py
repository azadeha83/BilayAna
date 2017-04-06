from setuptools import setup

setup(name='bilana',
      version='1.0',
      description='Analyse lipid bilayer systems',
#      url='',
      author='Fabian Keller',
      author_email='fabiankeller@wwu.de',
#      license='',
      packages=['bilana'],
      entry_points={
          'console_scripts': [
              'submit_energycalcs = bilana.command_line:submit_energycalcs',
              'initialize_system = bilana.command_line:initialize_system',
              'mend_energyruns = bilana.command_line:mend_energyruns',
              ]
          },
       #package_data={'': ['energy_recalculation.mdp']}
#      install_requires=['',],
#      scripts=['scripts/submit_energycalculations',],
      zip_safe=True)
