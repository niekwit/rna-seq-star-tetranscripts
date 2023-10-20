import pandas as pd
import click
import sys

def import_samples():
    try:
        csv = pd.read_csv("config/samples.csv")
        SAMPLES = csv["sample"]
        
        return SAMPLES

    except FileNotFoundError:
        click.secho("ERROR: config/samples.csv not found!", fg="red")
        sys.exit(1)