# -*- coding: utf-8 -*-
# Generated by Django 1.11 on 2016-07-25 16:26
from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('doseCalc', '0002_document'),
    ]

    operations = [
        migrations.DeleteModel(
            name='Document',
        ),
    ]
