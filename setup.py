from setuptools import setup

setup(name='src',
      version='1.0',
      description='Analyse lipid bilayer systems',
#      url='',
      author='Fabian Keller',
      author_email='fabiankeller@wwu.de',
#      license='',
      packages=['src'],
      entry_points={
          'console_scripts': [
              'submit_energycalcs = src.command_line:submit_energycalcs',
              'initialize_system = src.command_line:initialize_system',
              'mend_energyruns = src.command_line:mend_energyruns',
              'check_and_write = src.command_line:check_and_write',
              'write_eofscd = src.command_line:write_eofscd',
              'write_nofscd = src.command_line:write_nofscd'
              ]
          },
       #package_data={'': ['energy_recalculation.mdp']}
#      install_requires=['',],
#      scripts=['scripts/submit_energycalculations',],
      zip_safe=True)
