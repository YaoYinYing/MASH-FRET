---
layout: default
title: Getting started
nav_order: 1
---

<img src="assets/images/logos/logo-getting-started_400px.png" width="260" style="float:right; margin-left: 15px;"/>

# Getting started
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}


---

## System requirements

To run MASH-FRET, the software 
[MATLAB](https://fr.mathworks.com/products/matlab.html) and a list of specific toolboxes must be installed on your computer.
For more information about system requirements, please refer to the documentation section 
[System requirements](System_requirements.html).

**Note:** *MASH-FRET is tested to run under Matlab version R2011a and above*


---

## Installation
Clone or [download](https://github.com/RNA-FRETools/MASH-FRET/archive/master.zip) MASH-FRET into a directory of your choice.

```
git clone https://github.com/RNA-FRETools/MASH-FRET.git
```

Within MATLAB, add MASH-FRET to your search path by going to `Home → Set Path → Add with Subfolders`
{: .pb-3 }

Start MASH-FRET by typing `MASH` at the MATLAB command prompt.


---

## Interface

When you start MASH-FRET the **main user interface** and **action history** open simultaneously.

<a href="assets/images/gui/interface-default.png"><img src="assets/images/gui/interface-default.png" /></a>

The main user interface appears in its default layout and includes:
* a **tool bar** that allows the access the different <u>modules</u> (default: Video Processing) and to set the default <u>root folder</u> where files will be exported,
* a **main area** that contains the selected module,
* a **menu bar** that allows to modify the action history <u>view</u>, set <u>file overwriting</u> options and execute program <u>routines</u> (*advanced use*).

The action history lists the main events occurring in the program.
To keep track, logs are automatically written in a 
[daily log file](output-files/log-daily-logs.html).

**Note:** *Keep the action history or the command window of MATLAB visible when using MASH-FRET: information about running processes and failures are streaming in live.*


---

## Modules and tutorials

MASH-FRET's modules were originally created for the simulation, processing and analysis of video-based surface-immobilized single molecule FRET (smFRET) experiments.
Nonetheless, functionalities can be used for any type of video or trajectories requiring the same treatment.

The program offers five modules to work with:

- [Simulation](simulation.html)
- [Video processing](video-processing.html)
- [Trace processing](trace-processing.html)
- [Histogram analysis](histogram-analysis.html)
- [Transition analysis](transition-analysis.html)

Learn about module's individual functionalities in the respective documentation sections, or follow the 
[Tutorials](tutorials.html) to perform common tasks listed below:

* [Analyze experimental data](tutorials/analyze-data.html)
* [Validate analysis results](tutorials/validate-results.html)
* [Test algorithms](tutorials/test-algorithms.html)

