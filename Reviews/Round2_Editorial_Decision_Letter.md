From: Cerebral Cortex <onbehalfof@manuscriptcentral.com>  
Date: Mon, Nov 9, 2026 at 2:15 PM  
Subject: Cerebral Cortex - Manuscript Decision - CerCor-2022-00354R1  
To: <ssheldon@ualberta.ca>, <alona.fyshe@ualberta.ca>, <kyle.mathewson@ualberta.ca>  

09-Nov-2026  

**Manuscript ID:** CerCor-2022-00354R1  
**Title:** Connecting Covert Attention and Visual Perception to the Spatiotemporal Dynamics of Alpha Band Activity, Cross-Frequency Coupling (CFC), and Functional Connectivity using Multivariate Pattern Analysis (MVPA)  
**Authors:** Sarah S. Sheldon, Alona Fyshe, and Kyle E. Mathewson  

Dear Ms. Sheldon, Dr. Fyshe, and Dr. Mathewson,

Your revised manuscript (`CerCor-2022-00354R1`), along with your point-by-point response letter (`manuscript/response_to_reviewers.tex`) and open-source analysis/reproducibility pipeline (`OrientExoML`), has now been evaluated by the two original reviewers and discussed in a live editorial consultation between the reviewers and myself.

First, I would like to commend you for the substantial conceptual, structural, and visual improvements you made in response to the Round 1 reviews. Both reviewers were impressed by your expanded Discussion subsections addressing the study's novel contribution as a head-to-head single-trial multivariate benchmark across 10 oscillatory/CFC/connectivity estimators, single-trial estimator signal-to-noise ratio (SNR) differences, visuomotor and working-memory variance in the continuous $360^\circ$ orientation reproduction task, the visual field asymmetry in behavioral error magnitude, and the distinction between invasive ECoG (Voytek 2010), non-invasive MEG (Bonnefond and Jensen 2015), and scalp EEG high-gamma recordings. Furthermore, the restructured Introduction, formal academic tone, explicit methodological definitions (FOOOF $1/f$ vs. peak parameterization, IAF summary statistics, sideband-preserving gamma filter bandwidths, and HEOG step thresholds), high-resolution 300 DPI figures, and the newly added 3-panel summary comparison figure (**Figure 7**) have greatly strengthened the manuscript.

However, both reviewers independently conducted a detailed cross-check of the revised manuscript prose (`manuscript/main.tex` / `manuscript/generated_stats.tex`), the statistical JSON store (`manuscript/analysis_results.json`), the newly rendered figures (**Figures 2–7**), and the underlying Python reconstruction scripts (`manuscript/scripts/eeg_preproc_and_feats.py`, `manuscript/scripts/run_mvpa_and_figures.py`, and `manuscript/scripts/run_behavior_pipeline.py`), uncovering **critical internal contradictions between the manuscript text/captions and Figures 3–7**, as well as **several concrete code-versus-Methods discrepancies** that caused your Python pipeline to diverge from your original MATLAB pipeline.

Because the conceptual and pedagogical foundation of the paper is now strong, I am pleased to invite a **Minor-to-Moderate Technical Revision (Conditional Acceptance pending Text–Figure–Code Synchronization Verification)**. Before the manuscript can be formally accepted for publication in *Cerebral Cortex*, you must resolve the specific synchronization, methodological, and figure-formatting items itemized in the **Consolidated Editorial Action Checklist** below.

---

## Consolidated Editorial Action Checklist (Mandatory for Final Acceptance)

### 1. Resolve the "Split-Brain" Discrepancy Between Manuscript Text/Macros (`generated_stats.tex`) and Figures 3–7 (`analysis_results.json` Synchronization Bug)
When `manuscript/scripts/run_mvpa_and_figures.py` (lines 1636–1670) was executed to regenerate **Figures 3–6** and create **Figure 7**, it updated `analysis_results.json["mvpa_summary_metrics"]` and five keys in `analysis_results.json["mvpa"]["svm_results"]` (`alpha_pct_sig_range`, `alpha_pct_sig_min`, `alpha_left_sig_windows`, `alpha_right_sig_windows`, and `mvl_pct_sig_range`), **but left all remaining keys in `svm_results` and `svr_results` hardcoded to legacy numbers from an earlier run**. Consequently, the compiled manuscript text, figure captions, and response letter currently contradict the plotted figures in three major places:
1. **SVM Classification of Local CFC (Figure 3 vs. Text & Figure 7A–B):**
   - `main.tex` (lines 289–290) states that among the five local CFC metrics, *"only $\text{MVL}_{\text{\"Ozkurt}}$ yielded above-chance SVM classification... at $-875\text{ ms}$ pre-stimulus for left targets, and at $-785\text{ to }-770\text{ ms}$ and $-380\text{ ms}$ pre-stimulus for right targets"* and that *"None of the remaining CFC measures (AAC, MI, $\text{PAC}_{\text{GLM}}$, or $\text{PAC}_{\text{PLV}}$) classified covert attention above chance at any pre-stimulus time point."*
   - Yet in **Figure 3A** and **Figure 7A–B**, **MI** has green significance bars, inset topographies, and asterisks for both Left Targets ($-215\text{ ms}$, $2.0\%$) and Right Targets ($-500\text{ to }-485\text{ ms}$, $4.0\%$); **$\text{PAC}_{\text{PLV}}$** has green significance bars, inset topographies, and asterisks for both Left Targets ($-845\text{ ms}$, $2.0\%$) and Right Targets ($-425\text{ ms}$, $2.0\%$); **$\text{MVL}_{\text{\"Ozkurt}}$** has **4 separate green bars** on Left Targets ($8.0\%$) and **5 separate clusters (7 steps)** on Right Targets ($14.0\%$); and **Figure 3B** shows a green significance star at `Cz` for `PAC_GLM` (Left Targets) and only a single star near `CP1`/`C3` for `Alpha Amplitude` (Left Targets).
2. **SVM Classification of Functional Connectivity (Figure 4 vs. Text, Caption, & Figure 7A–B):**
   - `main.tex` (lines 302–312) and the Figure 4 caption state that **only wPLI** exceeded chance at a single window ($-462.5\text{ ms}$ on right-target trials) and that Correlation, PLV, and iCOH remained at chance throughout the pre-stimulus interval.
   - Yet in **Figure 4A** and **Figure 7A–B**, **wPLI has no green significance bar ($0.0\%$ on both Left and Right)**, whereas **Correlation (Left Targets, $-712.5\text{ ms}$, $8.3\%$)** and **PLV (Left Targets, $-412.5\text{ ms}$, $8.3\%$)** display green significance bars, 3D connectome insets, and blue asterisks in Figure 7A–B!
3. **SVR Prediction of Response Errors (Figures 5, 6, & 7C vs. Text & Captions):**
   - `main.tex` (lines 317–343) and the Figure 6 caption report single significant SVR windows for **Alpha Amplitude** ($-640\text{ ms}$, Non-Informative Right), **AAC** ($-200\text{ ms}$, Non-Informative Left), and **Correlation** ($-562.5\text{ ms}$, Non-Informative Left), while stating that *"None of the four phase–amplitude coupling metrics ($\text{MI}$, $\text{MVL}_{\text{\"Ozkurt}}$, $\text{PAC}_{\text{GLM}}$, or $\text{PAC}_{\text{PLV}}$) significantly predicted trial-by-trial response errors in any condition."*
   - Yet in **Figures 5A, 6A, and 7C**, Alpha Amplitude, AAC, and Correlation have **zero** significant bars or asterisks, and instead **$\text{PAC}_{\text{GLM}}$ (Informative Left, around $-800\text{ ms}$)** has a green significance bar, inset topography, and blue asterisk in Figure 7C!

**Required Action:** After fixing the Python signal-processing and statistical bugs detailed in Item 2 below, update `run_mvpa_and_figures.py` and `generate_latex_stats.py` so that **every single key** in `analysis_results.json["mvpa"]["svm_results"]`, `["svr_results"]`, and `["mvpa_summary_metrics"]` is populated dynamically from the exact same analysis run that renders **Figures 3–7**, and verify that the Results prose, Discussion, Figure 3–7 captions, and `response_to_reviewers.tex` match the plotted figures with 100% precision.

---

### 2. Fix Four Signal-Processing & Statistical Bugs in the Python Pipeline (`eeg_preproc_and_feats.py` and `run_mvpa_and_figures.py`)
Our reviewers identified the exact implementation bugs in the Python port that caused it to diverge from your original MATLAB pipeline (`repos/MVPA_CuedOrientTask_paper/Analysis_Code/`):
1. **Group-Level Permutation Baseline Variance Deflation (`run_mvpa_and_figures.py`, lines 310, 421, 1402–1404, 1444–1445):**
   - Your Methods text (line 254; `\statCVRepeatsNull{20}`, `\statCVTotalNull{200}`) and original MATLAB scripts (`permut_rep = 20`) specify **20 repetitions $\times$ 10 folds = 200 null evaluations** per participant, and DDTBOX (`Analyse_Decoding_Results_v2.m` line 223: `drawmode = 2`; `tmax_classifier_accuracies.m` lines 32–34) compares each participant's real cross-validated score against a **single drawn repetition** from that participant's permutation distribution so that the empirical and null series have matched estimation variance.
   - In `run_mvpa_and_figures.py`, only 10 null repetitions (100 evaluations) were run, and line 1402 (`all_perm_mean = np.mean(all_perm_reps, axis=(2, 3))`) averaged all 100 null shuffles per participant before computing `diff_scores = all_real - all_perm_mean`. Averaging 100 shuffles collapses the variance of `all_perm_mean` by a factor of $10$ ($\text{SD}$ shrinks by $\sqrt{10}$), artificially shrinking the standard error of `diff_scores`, inflating group $t$-statistics across all metrics, and causing the proliferation of false-positive single-window spikes in `MI`, `PAC_PLV`, `Correlation`, `PLV`, and `PAC_GLM` (as well as the nearly invisible permutation SEM error bars in Figures 3–6).
   - **Fix:** Restore 20 null repetitions (or match the exact DDTBOX `drawmode = 2` single-repetition draw / variance-matched null baseline from `tmax_classifier_accuracies.m`) so that the null baseline variance and group $t_{\max}$ tests are statistically valid.
2. **Biased Standard wPLI Instead of Vinck's Debiased Squared wPLI (`eeg_preproc_and_feats.py`, line 551):**
   - Both the manuscript text (Table 2, p. 11, p. 14) and your original MATLAB function (`conn_wPLI.m` line 64) specify the **debiased squared weighted phase-lag index** (Vinck et al. 2011):
     $$\text{dwPLI} = \frac{\left(\sum_t \text{Im}(X_t)\right)^2 - \sum_t \text{Im}(X_t)^2}{\left(\sum_t |\text{Im}(X_t)|\right)^2 - \sum_t \text{Im}(X_t)^2}$$
   - However, `eeg_preproc_and_feats.py` (line 551: `wpli_out[ij, pair_idx, tr] = abs(wpli_num) / wpli_den`) computes the **un-debiased standard wPLI**.
   - **Fix:** Implement the exact Vinck et al. (2011) debiased wPLI formula from `conn_wPLI.m` in `eeg_preproc_and_feats.py`.
3. **Missing Pre- and Post-EMCP $200\text{ ms}$ Baseline Subtractions (`eeg_preproc_and_feats.py`, lines 778–825):**
   - Both `main.tex` (line 182) and `Preproc_OrientWheel_Exo_ML.m` (lines 92 and 171) perform two $200\text{ ms}$ post-fixation baseline subtractions (`pop_rmbase`): one before the $\pm 1000\,\mu\text{V}$ threshold and Gratton ocular correction (EMCP), and a second one after EMCP before the $\pm 500\,\mu\text{V}$ threshold. `eeg_preproc_and_feats.py` omits both baseline subtractions, altering which epochs exceed the voltage thresholds due to DC offset.
   - **Fix:** Add both $200\text{ ms}$ post-fixation baseline subtractions in `eeg_preproc_and_feats.py` (and ensure the filter description in `main.tex` line 182 accurately reflects the preprocessing pipeline).
4. **Haufe Feature Weight Normalization (`/ 20.0` vs. `/ 100.0`), Inset Permutations (`2000` vs. `10000`), Behavioral Trial Filtering, and FOOOF Frequency Range:**
   - In `run_mvpa_and_figures.py`, fix the `/ 20.0` divisor in `ddtbox_svm_subject` (line 426) and `ddtbox_svr_subject` (line 584) to `/ 100.0` (matching the 10 repeats $\times$ 10 folds = 100 iterations) and verify Haufe pattern averaging across folds; change `n_perm=2000` to `n_perm=10000` on line 999 to match `main.tex` (line 259); ensure `main.tex` (line 184) and `run_behavior_pipeline.py` (lines 67–88) agree on whether behavioral ANOVAs include all completed behavioral trials or only EEG-artifact-free trials; and harmonize the FOOOF frequency range between `main.tex` (`3.9 to 48.8 Hz`, line 190) and `eeg_preproc_and_feats.py` (`2.0 to 50.0 Hz`, line 597).

---

### 3. Correct "Signed Response Errors" to "Absolute Response Errors" for SVR in the Manuscript Text
- In `main.tex` (page 16, line 252), the Methods text states: *"For support vector regression (SVR) ($\epsilon\text{-SVR}$, cost parameter $C = 0.1$), models were trained separately within each of the four Cue Type $\times$ Target Side conditions to predict continuous single-trial **signed response errors** from the pre-stimulus spatiotemporal EEG features."*
- However, both your original MATLAB code (`PrepData_DDTBOX.m`, lines 115–118: `SVR_labels{1,1} = squeeze(abs(errordeg...))`) and your Python code (`eeg_preproc_and_feats.py`, line 879: `c_svr_labels = np.abs(err_deg[mask])`) trained SVR on **absolute response error** ($|\text{Target} - \text{Response}| \in [0^\circ, 180^\circ]$).
- As both reviewers noted, predicting **absolute response error** (error magnitude / inverse perceptual precision) is the scientifically appropriate test of BSEM vs. OPR (since pre-stimulus oscillations prior to a randomly oriented target cannot predict clockwise vs. counterclockwise directional bias).
- **Fix:** Replace "signed response errors" with **"absolute response errors (orientation error magnitude)"** in `main.tex` (line 252) and clarify in the Results and Figure 5–7 captions that SVR predicted trial-by-trial absolute orientation response error.

---

### 4. Refine Table 1, Figure 2, Figures 3–6 Captions/Legends, and Minor Typos
1. **Table 1 (OPR Model Prediction 1, `main.tex` lines 123–126):**
   - Soften the OPR column under *1. Attentional Classification (SVM)* so that it does not claim the OPR model (Zazio et al. 2020) requires strictly *"At-chance classification of non-informative vs. informative cued trials by alpha amplitude alone."* Instead, state that under OPR, local alpha–gamma PAC or AAC should actively mediate and classify covert attentional states alongside alpha amplitude (and, critically under *2. Perceptual Error Prediction (SVR)*, that PAC or AAC must predict trial-by-trial perceptual error magnitude, whereas BSEM predicts that alpha amplitude classifies attentional state without predicting objective discrimination precision on non-informative trials).
2. **Figure 2 (`figure2_behavior.png` and `main.tex` lines 268–278):**
   - **Panel C (Staircase Threshold):** `figure2_behavior.png` plots three panels (`A. Absolute Orientation Error`, `B. Directional Bias`, and `C. Staircase Threshold`), but **Panel C** is never mentioned in the Figure 2 caption or the manuscript text. Please add a description of **(C)** to the Figure 2 caption and reference `(Figure 2C)` in the text where the staircase luminance distribution ($M = 54.0, \text{SD} = 22.0$) is reported.
   - **Significance Annotations on Panel A:** The Figure 2 caption references `(*p < 0.05, FDR-corrected tmax permutation tests; ns = not significant)`, but no `*` or `ns` brackets are drawn on Figure 2A. Please update `run_behavior_pipeline.py` to draw the `*` bracket over the Left Target bars and `ns` over the Right Target bars in Panel A.
3. **Figures 3–6 Captions, Legends, and Inset Positioning (`main.tex` lines 298, 311, 329, 342; `run_mvpa_and_figures.py`):**
   - **Trace Color Descriptions & Legend:** Update the captions of Figures 3, 4, 5, and 6 so they do not refer to the permuted-label chance traces as *"blue traces"* (in the rendered figures, permuted chance traces are plotted as **gray dashed traces**, while blue is used for the empirical Informative Cue trace in Figures 5 and 6). Add a compact legend to Panel A of Figures 3–6 so readers can immediately identify the empirical and permuted traces.
   - **Vertical vs. Horizontal Dashed Reference Line:** Either draw a vertical dashed line at $0\text{ ms}$ in Figures 3–6 or update the captions to refer to the horizontal dashed chance reference line ($50\%$ for SVM; $r_z = 0.0$ for SVR).
   - **Inset Topography Nose Clipping & Inset Count:** In `run_mvpa_and_figures.py` (line 1018), slightly lower the inset y-position (e.g., from `0.56` to `0.51`) so the nose tip on the circular scalp topographies in Figures 3A, 4A, and 5A is not clipped by the top subplot border, and ensure all significant time windows described in the text have corresponding visual representation.
4. **Citation Typo (`main.tex` line 250):**
   - Change `"itchie et al. 2019"` to `"Ritchie et al. 2019"`.

---

When submitting your revised manuscript (`CerCor-2022-00354R2`), please provide an updated point-by-point *Response to Reviewers* detailing how each item in the checklist above has been verified and resolved.

We look forward to receiving your final revised manuscript and are pleased to be able to move this comprehensive study toward publication in *Cerebral Cortex*.

Yours sincerely,

**Pasko Rakic, MD, PhD**  
Editor-in-Chief, *Cerebral Cortex*

---

# Appendix: Full Independent Round 2 Reviewer Reports

---

## Reviewer 1 — Round 2 Report

**Manuscript ID:** CerCor-2022-00354R1  
**Recommendation:** Minor / Major Revision (Pending Resolution of Text–Figure Statistical Synchronization and Methodological Discrepancies)

### General Assessment
In this revised manuscript (`CerCor-2022-00354R1`), the authors have put considerable effort into addressing the conceptual and interpretive concerns I raised in Round 1, as well as the presentation and methodological clarifications requested by Reviewer 2. In particular, the expanded Discussion sections—*"What the Current MVPA Comparison Adds to the Literature and Signal-to-Noise Considerations"*, *"Task Demands, Continuous Motor Reproduction, and Visual Field Asymmetry"*, and *"Implications for the OPR Model and Electrophysiological Scope of Scalp EEG"*—substantially improve the theoretical framing of the paper. The authors now clearly articulate their study's value as a head-to-head multivariate benchmark of first-order alpha amplitude against nine higher-order local cross-frequency coupling (CFC) and inter-electrode functional connectivity estimators within the same dataset and cross-validation framework. They also candidly discuss how single-trial estimator signal-to-noise ratio (SNR), continuous mouse-driven visuomotor reproduction variance, right-handedness and pooled staircase luminance, and within-analysis versus experiment-wise familywise Type-I error control shape the interpretation of their findings.

However, a rigorous audit of the revised manuscript text (`manuscript/main.tex`), the response letter (`manuscript/response_to_reviewers.tex`), the auto-generated statistical macros (`manuscript/generated_stats.tex` / `manuscript/analysis_results.json`), the newly reconstructed figures (`Figures 2–7`), and the underlying Python reconstruction scripts (`manuscript/scripts/eeg_preproc_and_feats.py` and `manuscript/scripts/run_mvpa_and_figures.py`) reveals **critical internal contradictions between the prose/macros and the newly generated figures**, as well as several **methodological discrepancies between what the Methods text states and what the analysis code actually executes**.

Specifically, when the authors re-ran the MVPA pipeline in Python (`run_mvpa_and_figures.py`) to regenerate high-resolution Figures 3–6 and create the new summary Figure 7, the set of metrics and time windows reaching statistical significance ($p < 0.05$, $t_{\max}$ corrected) changed relative to the original MATLAB analysis. While `analysis_results.json["mvpa_summary_metrics"]` and Figures 3–7 were updated with the new Python run's outputs, **half of `analysis_results.json["mvpa"]` (`svm_results` and `svr_results`) and the corresponding Results, Discussion, Figure Legends, and Response to Reviewers text were left hardcoded to the old single-window spikes!** Consequently, the manuscript text and figure captions currently describe green significance bars and asterisks that do not exist in Figures 3–7, while failing to mention the significance bars and asterisks that *are* plotted in Figures 3–7.

### Evaluation of Responses to Round 1 Comments
1. **Point 1 (Scalp EEG Limitations for Gamma/Connectivity & What the Study Adds Beyond Alpha Amplitude Decoding):** Satisfactory conceptually. The dedicated subsection (*"What the Current MVPA Comparison Adds to the Literature and Signal-to-Noise Considerations"*) and **Figure 7** provide a clear and defensible "take-home" message, and the discussion of invasive ECoG (Voytek 2010) versus MEG (Bonnefond and Jensen 2015) and scalp EEG attenuation/myogenic constraints above $50\text{ Hz}$ is well-reasoned.
2. **Point 2 (Signal-to-Noise Ratio [SNR] Confounds Across EEG Metrics):** Partially satisfactory, though qualitative. In `Table 1` (*Overview of predicted outcomes*), the authors state that the OPR model predicts *"At-chance classification of non-informative vs. informative cued trials by alpha amplitude alone."* The OPR model (Zazio et al. 2020) does not deny that alpha amplitude is modulated by spatial attention cues; rather, it proposes that alpha oscillations modulate sensory responsiveness via AAC (threshold) or PAC (gain). The authors should soften `Table 1` so as not to set up a "straw-man" prediction for the OPR model.
3. **Point 3 (Continuous $360^\circ$ Mouse-Rotation Task & Motor/Working Memory Variance):** Satisfactory conceptually in the Discussion, but undermined by a Methods vs. Code discrepancy: the manuscript text (line 252) claims SVR predicted *signed* response error ($[-180^\circ, +180^\circ]$), whereas the code (`eeg_preproc_and_feats.py` line 879) actually trained SVR on *absolute* response error ($|\text{error}| \in [0^\circ, 180^\circ]$).
4. **Point 4 (Multiple Comparison Correction Across 10 Metrics and Target Sides):** Conceptually sound in principle, but empirically contradicted by the new Python MVPA run (`run_mvpa_and_figures.py`), where unadjusted within-analysis $t_{\max}$ tests (inflated by the null-variance bug) produced spikes in 6 out of 10 metrics for Left Targets (`Alpha Amplitude`, `MVL_Özkurt`, `MI`, `PAC_PLV`, `FC: Correlation`, and `FC: PLV`), 4 out of 10 metrics for Right Targets, and an SVR spike in `PAC_GLM`—all ignored by the manuscript prose.
5. **Point 5 (Visual Field Asymmetry in Behavioral Performance, Figure 2):** Satisfactory. Both right-handed mouse control and the pooled 1-up/2-down detection staircase are clearly explained.

### Summary of Required Revisions (Reviewer 1)
1. **Fix the MVPA Null Comparison & Synchronize All Text, Macros, and Figures:** Correct the permutation null baseline in `run_mvpa_and_figures.py` so that it matches DDTBOX (`20 repeats` and proper null variance handling rather than averaging across 100 null shuffles), update all fields of `analysis_results.json["mvpa"]` and `analysis_results.json["mvpa_summary_metrics"]` from the same analysis run, and ensure that every significant window and asterisk shown in **Figures 3–7** matches the **Results text**, **Discussion**, **Figure Captions**, and **Response to Reviewers**.
2. **Reconcile Methods vs. Code Discrepancies:** Explicitly clarify that SVR predicted *absolute* response error (aligning `main.tex` line 252 and `eeg_preproc_and_feats.py` line 879), ensure behavioral trial rejection matches the Methods description (or state accurately how behavioral trials were filtered), ensure the two $200\text{ ms}$ baseline subtractions in `eeg_preproc_and_feats.py` match the EEG Preprocessing text, fix the `/ 20.0` vs. `/ 100.0` Haufe weight divisor and `n_perm=2000` vs. `10000` inset test, and harmonize the FOOOF frequency range (`3.9 to 48.8 Hz` vs. `2.0 to 50.0 Hz`).
3. **Harmonize Figure 2 and Figures 3–6 Captions:** Add `*` / `ns` annotations to Figure 2A and describe **Figure 2C** (staircase luminance distribution) in the Figure 2 caption and Results text; update the Figures 3–6 captions to match the gray dashed permutation traces and add a Panel A legend; and fix the typo in `Ritchie et al. 2019`.

---

## Reviewer 2 — Round 2 Report

**Manuscript ID:** CerCor-2022-00354R1  
**Recommendation:** Minor-to-Moderate Revision (Pending Resolution of Text–Figure–Code Synchronization Discrepancies)

### General Assessment
The authors have done a commendable job addressing the conceptual, stylistic, and pedagogical concerns I raised in Round 1:
1. **Clarity and Tone:** Colloquialisms ("the idea is", "As can be seen in...") have been removed, background theory is cleanly separated from the study's own hypotheses in the Introduction, and the neurophysiological discussion of invasive ECoG versus non-invasive MEG/EEG high-gamma recordings (Voytek 2010 vs. Bonnefond and Jensen 2015) is now nuanced and electrophysiologically sound.
2. **Methodological Transparency:** Explicit explanations are now provided for the post-mask fixation screen, self-paced breaks, signed vs. absolute orientation error definitions, the $25\,\mu\text{V}$ split-half sliding-window HEOG step criterion, FOOOF $1/f$ vs. peak parameterization, IAF summary statistics ($M = 10.6\text{ Hz}$, $\text{SD} = 0.5\text{ Hz}$), local within-electrode CFC scope, and the sideband-preservation rationale ($15\text{--}50\text{ Hz}$) for the gamma bandpass filter.
3. **Visual Overhaul and New Summary Figure:** Figures 2–6 have been regenerated at 300 DPI with readable typography and $200\text{ ms}$ tick spacing, and the new **Figure 7** provides the exact cross-metric summary comparison requested in Round 1.

However, cross-checking `main.tex`, `analysis_results.json`, **Figures 2–7**, and the Python scripts (`eeg_preproc_and_feats.py` and `run_mvpa_and_figures.py`) against the original MATLAB repository (`repos/MVPA_CuedOrientTask_paper/Analysis_Code/`) revealed **four concrete implementation bugs in the Python port** and **direct contradictions between the prose/captions and Figures 3–7**:
1. **Direct Contradictions Between Prose (`main.tex` / `generated_stats.tex`) and Rendered Figures 3–7:** `run_mvpa_and_figures.py` updated `mvpa_summary_metrics` (driving Figure 7) and rendered Figures 3–6 from the Python run, but left stale numbers in `analysis_results.json["mvpa"]["svm_results"]` and `["svr_results"]` for MVL, MI, $\text{PAC}_{\text{PLV}}$, wPLI, Correlation, PLV, and all SVR metrics.
2. **Four Signal-Processing & Statistical Bugs in the Python Port:**
   - **Biased Standard wPLI Instead of Debiased wPLI (`eeg_preproc_and_feats.py` line 551):** Computes $|\sum \text{Im}(X_t)| / \sum |\text{Im}(X_t)|$ instead of Vinck et al.'s (2011) debiased squared wPLI (`conn_wPLI.m` line 64).
   - **Group-Level Permutation Baseline Variance Deflation (`run_mvpa_and_figures.py` lines 310, 421, 1402–1404):** Computes 10 null repeats (100 evaluations) instead of 20 repeats (200 evaluations) and subtracts the mean across all 100 null evaluations (`all_perm_mean`) instead of drawing a single repetition per subject (`drawmode = 2` in `Analyse_Decoding_Results_v2.m` / `tmax_classifier_accuracies.m`), deflating baseline variance by $\sim 10\times$ and inflating group $t$-statistics.
   - **Missing Pre- and Post-EMCP Baseline Corrections (`eeg_preproc_and_feats.py` lines 778–825):** Omits the two $200\text{ ms}$ baseline subtractions (`pop_rmbase` in `Preproc_OrientWheel_Exo_ML.m` lines 92 & 171) prior to the $\pm 1000\,\mu\text{V}$ and $\pm 500\,\mu\text{V}$ rejection thresholds.
   - **FOOOF Frequency Range (`eeg_preproc_and_feats.py` lines 595–660 vs. `main.tex` p. 10):** Harmonize `3.9 to 48.8 Hz` vs. `2.0 to 50.0 Hz`.
3. **SVR Target Variable (Signed vs. Absolute Response Error):** Both `PrepData_DDTBOX.m` (line 115) and `eeg_preproc_and_feats.py` (line 879) train SVR on **absolute response error** (`abs(err_deg)`), whereas `main.tex` (line 252) erroneously states **"signed response errors"**. Correct the text to **"absolute response errors"**.
4. **Figure 2 & Figures 3–6 Formatting and Caption Fixes:** Describe **Figure 2C** (Staircase Threshold) in the text/caption and add `*` / `ns` brackets to **Figure 2A**; update the **Figures 3–6** captions to refer to **gray dashed traces** (not "blue traces") for the permuted-label null baseline and add a Panel A legend; prevent nose-tip clipping on the circular inset topographies in Figures 3A, 4A, and 5A; and fix `"itchie et al. 2019"` $\rightarrow$ `"Ritchie et al. 2019"`.
