# Metagenomics RD Swine

This repository contains the analysis pipeline and supporting scripts for a master's thesis project focused on the **Metagenomic characterization of swine nasal swab samples** using **Oxford Nanopore long-read sequencing**. 

## Overview
This project analyzes a subset of samples from a broader **Influenza A Virus (IAV) surveillance program** in swine. These particular samples were **negative for IAV by PCR** but came from animals that still showed **clinical signs of respiratory disease**. This raised the question of whether other viral agents might be involved.

The goal of the project is to identify and recover viral genomes from swine nasal samples using NGS metageomic data. The analysis pipeline includes **pre-processing, taxonomic classification and diversity analysis, viral genome recovery, and validation**, with steps tailored to maximize sensitivity while maintaining precision.

> **Before building the final pipeline**, a **benchmarking analysis** of multiple **assembly** and **classification tools** was performed using simulated and real datasets. This supported the selection of tools that **maximized recall and precision**, especially for low-abundance viral genomes in host-rich samples.
