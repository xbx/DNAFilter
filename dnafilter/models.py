from django.db import models

class SequenceBase(models.Model):
    db_name = models.CharField(max_length=50)
    db_path = models.CharField(max_length=300)
    db_version = models.CharField(max_length=15)
    def __unicode__(self):
        return self.db_name
