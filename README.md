DNAFilter
=========

Setup
-----

# Creates the database tables

django-admin.py syncdb
Creates the database tables for all apps in INSTALLED_APPS whose tables have not already been created.

Use this command when you've added new applications to your project and want to install them in the database. 
This includes any apps shipped with Django that might be in INSTALLED_APPS by default. 
When you start a new project, run this command to install the default apps.