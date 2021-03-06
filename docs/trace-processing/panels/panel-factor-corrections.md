---
layout: default
title: Factor corrections
parent: /trace-processing/panels.html
grand_parent: /trace-processing.html
nav_order: 6
---

<img src="../../assets/images/logos/logo-trace-processing_400px.png" width="170" style="float:right; margin-left: 15px;"/>

# Factor corrections
{: .no_toc }

Factor corrections is the fifth panel of module Trace processing.

Use this panel to apply cross-talk and gamma corrections.

<a class="plain" href="../../assets/images/gui/TP-panel-factors.png"><img src="../../assets/images/gui/TP-panel-factors.png" style="max-width: 292px;"/></a>

## Panel components
{: .no_toc .text-delta }

1. TOC
{:toc}


---

## Cross-talks settings

Use this interface to define the bleedthrough and direct excitation coefficients.

<a class="plain" href="../../assets/images/gui/TP-panel-factors-crosstalks.png"><img src="../../assets/images/gui/TP-panel-factors-crosstalks.png" style="max-width: 279px;"/></a>

The fluorescence of each emitter used in the experiment has the propensity to "bleed" trhough all unspecific detection channels. 
The ratio of fluorescence signal detected in the unspecific channel over the signel detected in the specific channel is called the bleedthrough coefficient.

To set the bleedthrough coefficient of a particular emitter in a particular unspecific channel, select the emitter in menu **(a)**, the destination channel in menu **(b)** and set the coefficient in **(c)**.

Beside, each emitter used in the experiment has the propensity to be excited by all unspecific lasers.
The ratio of fluorescence signal detected upon unspecific illumination over the signal detected upon specific illumination is called the direct excitation coefficient.

To set the direct excitation coefficient of a particular emitter by a particular unspecific laser, select the emitter in menu **(a)**, the unspecific laser in menu **(d)** and set the coefficient in **(e)**.

For more information about how to obtain cross-talk coefficients and how cross-talk correction is performed, see 
[Cross-talk correction](../workflow.html#cross-talk-correction) in Trace processing workflow.


---

## Gamma factor settings

Use this interface to define the gamma factors.

<a class="plain" href="../../assets/images/gui/TP-panel-factors-gamma.png"><img src="../../assets/images/gui/TP-panel-factors-gamma.png" style="max-width: 242px;"/></a>

Usually, the detection efficiency and quantum yield are different for the donor and the acceptor of a FRET pair. 
To put donor and acceptor intensities on the same scale and allow the conversion of apparent FRET to absolute distances, these differences must be corrected.
This is done by adjusting apparent FRET values with the gamma factor.

Gamma factors are specific to each FRET pair listed in **(a)** and can be set in three different ways:

- calculated from acceptor photobleaching
[<sup>1</sup>](#references) after selecting `From acceptor photobleaching` in menu **(b)**; see 
[Photobleaching-based calculation](#photobleaching-based-calculation) for more details
- set manually in **(c)** after selecting `Manual` in menu **(b)**
- loaded from an external ASCII file for all molecules after selecting `Manual` in menu **(b)** and pressing 
![Load](../../assets/images/gui/TP-but-load.png); see 
[.gam file](../../output-files/gam-gamma-factors.html) for more information about the file structure

For more information about how FRET values are corrected with gamma factors, see 
[Correct FRET values](../workflow.html#correct-fret-values) in Trace processing workflow.


### Photobleaching-based calculation
{: .no_toc }

Photobleaching-based gamma calculation requires intensity state trajectories; see 
[Data to discretize](panel-find-states.html#data-to-discretize) for more information about how to obtain state trajectories for intensity-time traces.

The gamma factor 
[*&#947;*<sub>*D*,*A*</sub>(*n*)](){: .math_var } correcting apparent FRET values for the donor-acceptor pair 
([*D*](){: .math_var },[*A*](){: .math_var }) on molecule 
[*n*](){: .math_var } can be calculated from intensity state trajectories 
[*state*<sub>*A*,em</sub><sup>*D*,ex</sup>](){: .math_var } and 
[*state*<sub>*D*,em</sub><sup>*D*,ex</sup>](){: .math_var }, before and after acceptor photobleaching such as:

{: .equation }
<img src="../../assets/images/equations/TP-eq-gamma-pb.gif" alt="\gamma_{D,A}(n) = \frac{state_{A,em}^{D,ex}(n,t_{0}-3)-state_{A,em}^{D,ex}(n,t_{0}+3)}{state_{D,em}^{D,ex}(n,t_{0}-3)-state_{D,em}^{D,ex}(n,t_{0}+3)}">

with 
[*t*<sub>0</sub>](){: .math_var } the detected photobleaching cutoff.

Acceptor photobleaching is detected in the acceptor intensity-time trace according to the Gamma factor options accessible by pressing 
![Opt.](../../assets/images/gui/TP-but-optp.png).

<a class="plain" href="../../assets/images/gui/TP-panel-factors-gamma-opt.png"><img src="../../assets/images/gui/TP-panel-factors-gamma-opt.png" style="max-width: 461px;"/></a>

Photobleaching detection settings are specific to each FRET pair selected in menu **(a)**.

Acceptor photobleaching is detected when the corresponding intensity-time trace drops below a certain intensity threshold defined in **(c)** providing a minimum cutoff value set in **(e)**.
To ensure detection at the very beginning of acceptor photobleaching, the detected cutoff position can be shifted downwards by a certain number of frames set in **(d)**.

When photobleaching was successfully detected, the icon 
![valid](../../assets/images/gui/TP-panel-factors-gamma-opt-valid.png "Valid") appears, otherwise the icon 
![invalid](../../assets/images/gui/TP-panel-factors-gamma-opt-invalid.png "Invalid") is shown to inform that gamma factors can not be computed.

Photobleaching cutoff is displayed in **(b)** in seconds or frame according to time-axis units defined in 
[Time axis](panel-plot.html#time-axis) and can be shown in the top axes with a red cursor by activating the option in **(f)**.

To save options and compute the FRET-pair-specific gamma factors, press 
![compute &#947;](../../assets/images/gui/TP-but-compute-gamma.png "compute &#947;").

Gamma factors can also be loaded from an external file for all molecules by pressing 
![load gamma](../../assets/images/gui/TP-but-load-gamma.png "load gamma"); see 
[.gam file](../../output-files/gam-gamma-factors.html) for more information about the file structure.


### References
{: .no_toc }

1. J.J. McCann, U.B. Choi, L. Zheng, K. Weninger, M.E. Bowen, *Optimizing Methods to Recover Absolute FRET Efficiency from Immobilized Single Molecules*, *Biophys J.* **2010**, DOI: [10.1016/j.bpj.2010.04.063](https://dx.doi.org/10.1016%2Fj.bpj.2010.04.063).


---

## Apply settings to all molecules

Press 
![all](../../assets/images/gui/TP-but-all.png "all") to apply emitter- and illumination-specific 
[Cross-talks settings](#cross-talks-settings), as well as the FRET pair-specific 
[Gamma factor settings](#gamma-factor-settings) to all molecules.

Corrections are applied to other molecules only when the corresponding data is processed, *e. g.* when pressing 
![UPDATE ALL](../../assets/images/gui/TP-but-update-all.png "UPDATE ALL"); see 
[Process all molecules data](panel-sample-management.html#process-all-molecules-data) for more information.

**Warning:** *Applying settings to all molecules will overwrite gamma factors imported from file with the factors defined for the current molecule. 
To recover the gamma factors from file, the gamma factor file must be re-imported; see 
[Gamma factor settings](#gamma-factor-settings) for help.*

