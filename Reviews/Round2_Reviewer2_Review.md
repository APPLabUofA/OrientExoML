# Cerebral Cortex — Round 2 Peer Review (Reviewer 2)

**Manuscript ID:** CerCor-2022-00354R1  
**Title:** Connecting Covert Attention and Visual Perception to the Spatiotemporal Dynamics of Alpha Band Activity, Cross-Frequency Coupling (CFC), and Functional Connectivity using Multivariate Pattern Analysis (MVPA)  
**Authors:** Sarah S. Sheldon, Alona Fyshe, and Kyle E. Mathewson  
**Recommendation:** Minor-to-Moderate Revision (Pending Resolution of Text–Figure–Code Synchronization Discrepancies)

---

## General Assessment

The authors have done a commendable job addressing the conceptual, stylistic, and pedagogical concerns I raised in Round 1. In particular:

1. **Clarity and Tone:** The writing throughout the Abstract, Introduction, Methods, Results, and Discussion is much more formal, precise, and readable. Colloquialisms ("the idea is", "As can be seen in...") have been removed, background theory is now cleanly separated from the study's own hypotheses and analytic approach in the Introduction, and the neurophysiological discussion of invasive ECoG versus non-invasive MEG/EEG high-gamma recordings (Voytek et al., 2010 vs. Bonnefond & Jensen, 2015) is now nuanced and electrophysiologically sound.
2. **Methodological Transparency:** The Methods section now provides explicit explanations of the post-mask fixation screen, self-paced block breaks, mathematical definitions of signed and absolute orientation error, the $25\,\mu\text{V}$ split-half sliding-window HEOG step criterion for excluding overt saccades, FOOOF parameterization of $1/f$-like aperiodic background activity versus Gaussian oscillatory peaks, descriptive statistics and filtering usage of Individual Alpha Frequency (IAF), local within-electrode CFC scope, and the sideband-preservation rationale (and upper-beta caveat) for the $15\text{--}50\text{ Hz}$ high-frequency bandpass filter.
3. **Visual Overhaul and New Summary Figure:** Figures 2–6 have been completely regenerated at 300 DPI with readable typography, $200\text{ ms}$ x-axis tick spacing, concise panel titles, and colorblind-safe palettes. The newly added **Figure 7** provides the exact cross-metric summary comparison of peak decoding accuracy, proportion of significant pre-stimulus time windows, and peak SVR correlation that I requested in Round 1.

However, a rigorous cross-check of the revised manuscript text (`main.tex` / `generated_stats.tex`), the statistical summary file (`analysis_results.json`), the newly rendered figures (`figure2_behavior.png` through `figure7_mvpa_summary.png`), and the accompanying Python analysis/plotting scripts (`eeg_preproc_and_feats.py` and `run_mvpa_and_figures.py`) revealed **several critical discrepancies between the manuscript prose/captions and the rendered figures/code**. Specifically, when the authors re-ran the pipeline in Python to generate Figures 3–7, `run_mvpa_and_figures.py` updated the summary arrays plotted in Figures 3–7 (`mvpa_summary_metrics`, `alpha_left_sig_windows`, `alpha_right_sig_windows`, and `mvl_pct_sig_range`), **but left stale numbers from an earlier run in `analysis_results.json` for all other SVM and SVR metrics**, which were then macro-expanded into the Results and Discussion prose. Furthermore, four subtle signal-processing and statistical differences between the Python reconstruction (`eeg_preproc_and_feats.py` / `run_mvpa_and_figures.py`) and the original MATLAB pipeline (`repos/MVPA_CuedOrientTask_paper/Analysis_Code/`) explain why the plotted significance masks in Figures 3–7 shifted from the numbers reported in the text.

Below, I first document the status of all Round 1 comments and then detail the exact text–figure–code discrepancies that must be resolved before publication.

---

## Part 1: Point-by-Point Evaluation of Round 1 Comments

### Abstract
- **Listing of three alpha functions too elaborate:** *Resolved.* Condensed into a single concise opening sentence.
- **Spell out CFC, SVM, and SVR at first use:** *Resolved.* All three acronyms are now defined at first mention in the Abstract.
- **"pre-target" $\rightarrow$ "pre-stimulus":** *Resolved.* Standardized to "pre-stimulus" throughout the Abstract and main text.
- **"metrics have activity" $\rightarrow$ "time courses" / "spatiotemporal patterns":** *Resolved.* Rephrased accurately.
- **"alpha amplitude alone is most sensitive":** *Resolved.* Revised to *"Among the tested metrics, pre-stimulus alpha amplitude was most sensitive to the deployment of covert spatial attention."*

### Introduction
- **Separate literature review from own approach/hypotheses (pages 3–6):** *Resolved.* The first four paragraphs and the unified theoretical comparison paragraph now review the BSEM, OPR, and functional connectivity literature, while the final paragraph of the Introduction (together with Table 1) cleanly introduces the present study's MVPA design and competing predictions.
- **Erroneous question mark on page 6 ("...information? (Jensen and Mazaheri 2010)"):** *Resolved.* Fixed.
- **Avoid colloquial "the idea is" (page 7):** *Resolved.* Replaced with formal academic phrasing.

### Materials and Methods
- **"identical to the fixation" (page 9):** *Resolved.* Explicitly clarifies that during the $516.6\text{ ms}$ post-mask delay, the display returned to the initial fixation configuration (central white fixation dot and two empty dark gray placeholder circles).
- **"rest at their own pace" (page 9):** *Resolved.* Clarifies that participants could pause between blocks of 48 trials for as long as they wished before pressing the spacebar to continue.
- **Signed vs. absolute response error & $F_{\max}$ ANOVA formulation (page 10):** *Resolved.* Signed error ($\theta_{\text{response}} - \theta_{\text{target}} \in [-90^\circ, +90^\circ)$) and absolute error ($|\theta_{\text{response}} - \theta_{\text{target}}| \in [0^\circ, 90^\circ]$) are now explicitly defined, and the permutation $F_{\max}$ description is clear.
- **How EOG rejection implies failure to use covert attention (page 11):** *Resolved.* The text now explains the split-half sliding-window step threshold ($25\,\mu\text{V}$ over $200\text{ ms}$ windows, $10\text{ ms}$ step, corresponding to $>1.5^\circ\text{--}2^\circ$ horizontal rotation) and why losing $>30\%$ of trials to this criterion indicates systematic overt eye movements toward the lateralized placeholders rather than covert attention.
- **FOOOF $1/f$ vs. peak components, IAF reporting, and IAF usage (page 12):** *Partially Resolved (see Major Issue 5 below for a minor code/text parameter discrepancy).* The prose now clearly explains that FOOOF models the power spectrum as the sum of a broadband aperiodic $1/f$-like decay ($L(F) = b - \log(F^\chi)$) and $N$ narrowband Gaussian peaks, reports the grand-average IAF ($M = 10.6\text{ Hz}$, $\text{SD} = 0.5\text{ Hz}$, range $= 9.6\text{--}11.4\text{ Hz}$), and explains how the condition- and electrode-specific median IAF centered the $\pm 1\text{ Hz}$ bandpass filters for functional connectivity.
- **Why phase–phase CFC was excluded & local vs. inter-areal CFC (page 12):** *Resolved.* The authors clearly justify focusing on local within-electrode AAC and PAC because the OPR hypothesis specifically predicts that local alpha phase/power gates local gamma duty cycles within task-relevant visual columns, while noting inter-areal CFC and phase–phase coupling as future extensions.
- **High-frequency gamma filtering (`Af1`, `Af2`, $15\text{--}60\text{ Hz}$ vs. $30\text{--}45\text{ Hz}$, and beta contamination) (page 13):** *Resolved.* The text now explicitly states that a single wide bandpass filter was applied (`Af1 = round(30 - (37.5 / 2.5)) = 15 Hz` to `Af2 = round(45 + (37.5 / 2.5)) = 60 Hz`, bounded above by the $50\text{ Hz}$ preprocessing low-pass filter, yielding an effective $15\text{--}50\text{ Hz}$ passband) in order to preserve the $\pm\text{alpha}$ sidebands required to detect alpha-rate amplitude modulation of the gamma envelope (citing Aru et al., 2015 and Seymour et al., 2017), and candidly acknowledges potential upper-beta ($15\text{--}30\text{ Hz}$) leakage.

### Results & Discussion
- **Avoid "As can be seen in..." (pages 17–21):** *Resolved.* Removed throughout.
- **Match text order of absolute vs. signed response error to Figure 2 (page 17):** *Resolved.* Absolute error (Figure 2A) is now presented before signed error (Figure 2B).
- **"during the 1000 ms before target onset" & "which predicts that changes" (page 21):** *Resolved.* Both sentences have been revised clearly.
- **Voytek et al. (2010) ECoG vs. Bonnefond & Jensen (2015) MEG vs. scalp EEG passage (page 23):** *Resolved.* The rewritten paragraph accurately explains that Voytek et al. (2010) observed a shift from theta–high-gamma ($80\text{--}150\text{ Hz}$) to alpha–high-gamma PAC in posterior subdural ECoG electrodes during visual tasks, that Bonnefond & Jensen (2015) demonstrated non-invasive MEG alpha–high-gamma ($80\text{--}120\text{ Hz}$) PAC outside the skull, and that the present study restricted scalp EEG analysis to $30\text{--}45\text{ Hz}$ because of cranial tissue attenuation, larger spatial integration volumes, and broadband cranial/ocular muscle artifacts above $50\text{ Hz}$ in scalp EEG.

---

## Part 2: Major Technical & Text–Figure–Code Discrepancies to Resolve in Round 2

### Major Issue 1: Direct Contradictions Between the Manuscript Prose (`main.tex` / `generated_stats.tex`) and Rendered Figures 3–7
When inspecting `run_mvpa_and_figures.py` (lines 1612–1660) alongside `analysis_results.json` and `generate_latex_stats.py`, it becomes clear why the text and figures disagree: `run_mvpa_and_figures.py` overwrote `mvpa_summary_metrics` (which drives Figure 7) and three specific keys (`alpha_left_sig_windows`, `alpha_right_sig_windows`, and `mvl_pct_sig_range`), **but did not update the remaining keys in `analysis_results.json["mvpa"]["svm_results"]` and `["svr_results"]`**. As a result, Figures 3–7 display the output of the new Python run (`data/cache/group_mvpa_results.npz`), whereas the manuscript prose and captions still cite legacy numbers from an earlier run:

1. **SVM Decoding of Local Alpha–Gamma CFC (Results p. 13 vs. Figure 3A, 3B & Figure 7A–B):**
   - **What the text claims (`main.tex` lines 278–285):** *"Among the five local alpha–gamma CFC metrics, only $\text{MVL}_{\text{\"Ozkurt}}$ yielded above-chance SVM classification... Specifically, $\text{MVL}_{\text{\"Ozkurt}}$ exceeded the permuted-label null distribution at $-875\text{ ms}$ for left targets, and at $-785\text{ to }-770\text{ ms}$ and $-380\text{ ms}$ for right targets ($8.0\text{--}14.0\%$ of pre-stimulus windows; Figure 3A; Figure 7A–B). None of the remaining CFC measures (AAC, MI, $\text{PAC}_{\text{GLM}}$, or $\text{PAC}_{\text{PLV}}$) classified covert attention above chance at any pre-stimulus time point."*
   - **What Figure 3A, Figure 3B, and Figure 7A–B actually show:**
     - **MI** *is* statistically significant (green shaded bars + inset topographies in Figure 3A; asterisks in Figure 7A; $2.0\%$ Left / $4.0\%$ Right bars in Figure 7B) at $-215\text{ ms}$ (step 44) for Left Targets and $-500\text{ to }-485\text{ ms}$ (steps 25–26) for Right Targets.
     - **$\text{PAC}_{\text{PLV}}$** *is* statistically significant (green shaded bars + inset topographies in Figure 3A; asterisks in Figure 7A; $2.0\%$ Left / $2.0\%$ Right bars in Figure 7B) at $-845\text{ ms}$ (step 2) for Left Targets and $-425\text{ ms}$ (step 30) for Right Targets.
     - **$\text{MVL}_{\text{\"Ozkurt}}$** in Figure 3A shows **4 separate significant windows** for Left Targets ($-875, -695, -560, \text{ and }-395\text{ ms}$, which is $4/50 = 8.0\%$) and **7 significant steps across 5 clusters** for Right Targets ($-800, -725, -470, -440, \text{ and }-380\text{ to }-350\text{ ms}$, which is $7/50 = 14.0\%$)—not a single window at $-875\text{ ms}$ and two windows at $-785\text{ to }-770\text{ ms}$ and $-380\text{ ms}$. (Notice that the text printed the new percentage range `$8.0--14.0%$` alongside the old window timestamps!)
     - In **Figure 3B**, the `PAC_GLM` Left Target topography displays a green significance star on central electrode `Cz`, even though the prose only reports significant epoch-averaged feature weights for Alpha Amplitude.

2. **SVM Decoding of Inter-Electrode Alpha Functional Connectivity (Results p. 14 & Figure 4 Caption vs. Figure 4A & Figure 7A–B):**
   - **What the text and caption claim (`main.tex` lines 298–309):** *"Among the four inter-electrode alpha functional connectivity metrics, only debiased weighted phase-lag index (wPLI) exceeded permuted-label chance at a single isolated time window ($-462.5\text{ ms}$ pre-stimulus for right targets; Figure 4A; Figure 7A–B)... Amplitude envelope correlation, PLV, and iCOH remained at chance throughout the entire pre-stimulus interval."* And the Figure 4 caption states: *"The green shaded bar indicates the single significant time window for wPLI on right-target trials ($-462.5\text{ ms}$)..."*
   - **What Figure 4A and Figure 7A–B actually show:**
     - **wPLI** has **no significant windows at all** ($0\%$ on both Left and Right targets; no green bar in Figure 4A, no asterisk in Figure 7A, $0\%$ bar in Figure 7B).
     - Instead, **Correlation** (Left Targets at $-712.5\text{ ms}$, $8.3\%$ of windows) and **PLV** (Left Targets at $-412.5\text{ ms}$, $8.3\%$ of windows) show green shaded significant bars and inset 3D connectomes in Figure 4A, as well as asterisks and $8.3\%$ bars in Figure 7A–B!

3. **SVR Prediction of Single-Trial Response Errors (Results p. 14–15 & Figure 6 Caption vs. Figures 5, 6, & 7C):**
   - **What the text and captions claim (`main.tex` lines 313–345):**
     - Alpha amplitude significantly predicted response error at $-640\text{ ms}$ on non-informative right-target trials (with significant occipital/parietal feature weights in Figure 5B).
     - AAC significantly predicted response error at $-200\text{ ms}$ on non-informative left-target trials (with significant fronto-central feature weights at Fz, FCz, FC1).
     - *"None of the four phase–amplitude coupling metrics (MI, $\text{MVL}_{\text{\"Ozkurt}}$, $\text{PAC}_{\text{GLM}}$, or $\text{PAC}_{\text{PLV}}$) significantly predicted trial-by-trial response errors in any condition."*
     - Functional connectivity **Correlation** significantly predicted response error at $-562.5\text{ ms}$ on non-informative left-target trials (and the Figure 6 caption explicitly directs the reader to *"The green shaded bar marks the single significant prediction window for amplitude envelope correlation ($-562.5\text{ ms}$ on non-informative left-target trials)"*).
   - **What Figures 5, 6, and 7C actually show:**
     - In **Figure 5A–B**, neither Alpha Amplitude nor AAC has any green significant bar or significant feature weight star in any condition! Instead, **$\text{PAC}_{\text{GLM}}$ on Informative Left trials** has a green significant bar and inset topography at $-800\text{ ms}$ (step 5), and **Figure 7C** places a blue significance asterisk on `CFC: PAC (GLM)` (Informative Cue)!
     - In **Figure 6A**, there is **no green shaded bar anywhere** in any of the eight panels (and no asterisk on any FC metric in Figure 7C), directly contradicting both the Results text and the Figure 6 caption!

---

### Major Issue 2: Four Signal-Processing & Statistical Bugs in the Python Reconstruction (`eeg_preproc_and_feats.py` and `run_mvpa_and_figures.py`)
Inspecting the original MATLAB codebase (`repos/MVPA_CuedOrientTask_paper/Analysis_Code/`) alongside the new Python scripts (`manuscript/scripts/eeg_preproc_and_feats.py` and `run_mvpa_and_figures.py`) reveals four concrete implementation bugs in the Python port that explain why the results diverged from the original MATLAB analysis:

1. **Biased Standard wPLI Instead of Debiased wPLI (`eeg_preproc_and_feats.py` line 551):**
   - The manuscript (Table 2, p. 11, p. 14) and the original MATLAB function (`conn_wPLI.m` line 64, called by `PrepConnData_DDTBOX.m` line 21) compute the **debiased squared weighted phase-lag index** (Vinck et al., 2011):
     $$\text{dwPLI} = \frac{\left(\sum_t \text{Im}(X_t)\right)^2 - \sum_t \text{Im}(X_t)^2}{\left(\sum_t |\text{Im}(X_t)|\right)^2 - \sum_t \text{Im}(X_t)^2}$$
   - However, `eeg_preproc_and_feats.py` (line 551: `wpli_out[ij, pair_idx, tr] = abs(wpli_num) / wpli_den`) computes the **un-debiased standard wPLI** ($|\sum \text{Im}(X_t)| / \sum |\text{Im}(X_t)|$). This directly explains why the wPLI decoding trace changed between the original MATLAB run and the Python figure regeneration.

2. **Group-Level Permutation Baseline Variance Deflation (`run_mvpa_and_figures.py` lines 310, 421, 1402–1404):**
   - The manuscript (p. 12) states that empirical chance baselines were generated using **20 repetitions** of 10-fold cross-validation with shuffled labels (200 null evaluations per participant, matching `permut_rep = 20` in `Run_ConnDecoding_v1.m` and `Run_SVR_v2.m`), and DDTBOX's group-level $t_{\max}$ test (`Analyse_Decoding_Results_v2.m` line 223: `drawmode = 2`; `tmax_classifier_accuracies.m` lines 32–34) compares each participant's real cross-validated accuracy against **one randomly drawn repetition** from that participant's null distribution (`all_subj_perm_acc_reps_draw`) so that the real and null scores have matched estimation variance.
   - In `run_mvpa_and_figures.py` (lines 310 & 421), only **10 null repetitions** (1 per CV repeat = 100 null evaluations) are computed, and lines 1402–1404 (`diff_scores = all_real - all_perm_mean`) subtract the **mean across all 100 null evaluations** (`all_perm_mean`). Because averaging 100 null evaluations reduces the variance of the baseline term by $\sim 10\times$, the standard error of `diff_scores` across participants is artificially deflated, inflating the group $t$-statistics and producing the isolated single-step false-positive spikes seen in Figures 3A, 4A, and 5A (`MI`, `PAC_PLV`, `Correlation`, `PLV`, `PAC_GLM`).

3. **Missing Pre- and Post-EMCP Baseline Corrections (`eeg_preproc_and_feats.py` lines 778–825):**
   - The Methods text (p. 9) and the original MATLAB preprocessing script (`Preproc_OrientWheel_Exo_ML.m` lines 92 & 171) perform **two** $200\text{ ms}$ baseline subtractions (`pop_rmbase`): one prior to Gratton ocular correction (EMCP) and a second one after EMCP prior to the $\pm 500\,\mu\text{V}$ second-stage voltage threshold rejection.
   - `eeg_preproc_and_feats.py` (lines 778–825) omits both baseline subtractions around `gratton_emcp_numpy`. Because raw un-baselined EEG epochs have DC offsets, applying the $\pm 1000\,\mu\text{V}$ and $\pm 500\,\mu\text{V}$ peak thresholds without baseline subtraction alters which trials are rejected per participant.

4. **FOOOF Frequency Range and Estimator (`eeg_preproc_and_feats.py` lines 595–660 vs. `main.tex` p. 10):**
   - The Methods text (p. 10) states that FOOOF (`v1.0.0`) was fit across a **$3.9\text{ to }48.8\text{ Hz}$** frequency range (`max_n_peaks = 6`, `min_peak_height = 0.0`, `peak_threshold = 2.0`, `peak_width_limits = [2.0, 12.0] Hz`, `aperiodic_mode = 'fixed'`).
   - In the repository's original `fooof_analysis.py` (line 82), `fg.fit(freqs, psds, [2, 50])` was called with `[2, 50] Hz` (at $\Delta f = 1000/512 = 1.953125\text{ Hz}$, bins from $1.95\text{ to }48.83\text{ Hz}$), and `eeg_preproc_and_feats.py` (lines 595–660) replaces `fooof.FOOOFGroup` with a custom Numba least-squares heuristic over `2.0 <= freqs[k] <= 50.0`. The authors should ensure the exact frequency range reported in `main.tex` matches the code and use the actual `fooof` package (or verify that the custom Numba routine reproduces identical IAF estimates).

---

### Major Issue 3: Contradiction Between Manuscript Text and Code on SVR Target Variable (Signed vs. Absolute Response Error)
- **What `main.tex` states (p. 12, Multivariate Pattern Analysis, lines 231–233):** *"Second, to test whether pre-stimulus neural dynamics predicted behavioral precision on a trial-by-trial basis, linear support vector regression (SVR; $c = 0.1$) models were trained separately within each of the four Cue Type $\times$ Target Side conditions to predict continuous single-trial **signed response errors** from the pre-stimulus spatiotemporal EEG features."*
- **What both the original MATLAB code and Python code actually compute:**
  - In `repos/MVPA_CuedOrientTask_paper/Analysis_Code/PrepData_DDTBOX.m` (lines 115–118):
    ```matlab
    SVR_labels{1,1} = squeeze(abs(errordeg{1,isubj}(1,1,nonzeros(trials_used{1,1}))));
    ```
  - And in `manuscript/scripts/eeg_preproc_and_feats.py` (line 879):
    ```python
    c_svr_labels = np.abs(err_deg[mask])
    ```
- **Why this matters:** Both codebases train SVR to predict **absolute response error** ($|\theta_{\text{response}} - \theta_{\text{target}}|$), *not* signed response error! Moreover, from a theoretical standpoint, predicting **absolute orientation error** (the magnitude of perceptual deviation, i.e., inverse precision) is the scientifically appropriate test of whether pre-stimulus alpha/CFC/connectivity fluctuations modulate visual perceptual quality. Predicting signed clockwise vs. counterclockwise error from cue-locked activity before a randomly oriented Gabor appears would have no coherent physiological interpretation. Please correct "signed response errors" to **"absolute response errors"** in the Methods (p. 12) and wherever else applicable.

---

### Major Issue 4: Figure 2 and Figures 3–6 Caption and Visual Formatting Mismatches
While the visual quality of the figures is vastly improved over Round 1, several discrepancies between the figure images and their captions/text remain:

1. **Figure 2 (`figure2_behavior.png` vs. `main.tex` lines 253–273):**
   - **Unreferenced Panel C:** `figure2_behavior.png` contains **three panels**: `A. Absolute Orientation Error`, `B. Directional Bias`, and `C. Staircase Threshold` (plotting Contrast in % Color Value for Left vs. Right targets). However, neither the Results text (pp. 12–13) nor the Figure 2 caption (`main.tex` lines 269–271) mentions **Panel C** at all! Either describe Panel C in the Results text and Figure 2 caption, or remove Panel C if only Panels A and B are intended.
   - **Missing Significance Markers (`*` and `ns`) in Panel A:** The Figure 2 caption states *"($*p < 0.05$, FDR-corrected $t_{\max}$ permutation tests; $\text{ns} = \text{not significant}$)"*, and in Round 1 I specifically noted how informative Figure 2 is when the behavioral cueing effect is visually indicated. In the newly rendered `figure2_behavior.png`, there are **no `*` or `ns` brackets/annotations** drawn on Panel A (or Panels B–C). Please either add the `*` / `ns` brackets to Figure 2A or update the caption accordingly.

2. **Figures 3–6 (`figure3_svm_cfc.png` through `figure6_svr_fc.png` vs. Captions in `main.tex` lines 289–340):**
   - **"Blue traces" vs. Gray traces:** All four captions (Figures 3, 4, 5, and 6) state *"blue traces denote permuted-label chance accuracy"* (or *"blue traces show permuted-label null $r_z$"*). In the actual rendered figures, the permuted-label null traces are plotted as **gray dashed lines with triangles (`#7F7F7F`)**, while the real Informative Cue traces in Figures 5 and 6 are plotted in blue (`#0072B2`)! Describing the null traces as "blue" in the captions is therefore particularly confusing in Figures 5 and 6, where blue actually represents the Informative Cue condition. Please update all four captions to state **"gray dashed traces denote permuted-label chance..."** and add a small legend to Panel A of Figures 3–6 so the reader can immediately distinguish Real vs. Permuted (in Figures 3–4) and Informative vs. Non-Informative vs. Permuted (in Figures 5–6).
   - **"Vertical dashed line marks target onset (0 ms)":** All four captions state *"The vertical dashed line marks target onset (0 ms)"*, but no vertical dashed line is drawn at $0\text{ ms}$ in Figures 3–6 (only a horizontal dotted line at chance = $50\%$ or $0.0$).
   - **Inset Topography Clipping and Cap of 2 Insets per Panel:** In Figures 3A, 4A, and 5A, the top tip of the nose on the circular inset topographies (`inset_axes([..., 0.56, 0.26, 0.42])`) is slightly clipped by the upper axis border. In addition, in Figure 3A (`MVL_Özkurt` Left Targets), there are 4 green significant windows, but only the first 2 have inset topographies because `run_mvpa_and_figures.py` (line 1018) caps insets at `if run_idx < 2:`. Lowering the y-position of the inset slightly (e.g., from `0.56` to `0.52`) will prevent nose clipping.

---

## Part 3: Minor Textual Corrections

1. **Page 12, Line 250 (`main.tex`):** Typo in citation: `"itchie et al. 2019"` $\rightarrow$ `"Ritchie et al. 2019"`.
2. **Table 1 & Discussion Terminology:** In the Discussion (p. 17), when discussing the isolated/transient effects in the secondary CFC and functional connectivity metrics, make sure the discussion matches whichever final synchronized results are shown in Figures 3–7 once the Python/MATLAB pipeline discrepancy (Major Issues 1 & 2) is resolved.
