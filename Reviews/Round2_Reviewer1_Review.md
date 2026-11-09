# Round 2 Peer Review — Reviewer 1

**Manuscript ID:** CerCor-2022-00354R1  
**Title:** Connecting Covert Attention and Visual Perception to the Spatiotemporal Dynamics of Alpha Band Activity, Cross-Frequency Coupling (CFC), and Functional Connectivity using Multivariate Pattern Analysis (MVPA)  
**Authors:** Sarah S. Sheldon, Alona Fyshe, and Kyle E. Mathewson  
**Journal:** *Cerebral Cortex*  
**Recommendation:** Minor / Major Revision (Pending Resolution of Text–Figure Statistical Synchronization and Methodological Discrepancies)

---

## General Assessment

In this revised manuscript (`CerCor-2022-00354R1`), the authors have put considerable effort into addressing the conceptual and interpretive concerns I raised in Round 1, as well as the presentation and methodological clarifications requested by Reviewer 2. In particular, the expanded Discussion sections—*"What the Current MVPA Comparison Adds to the Literature and Signal-to-Noise Considerations"*, *"Task Demands, Continuous Motor Reproduction, and Visual Field Asymmetry"*, and *"Implications for the OPR Model and Electrophysiological Scope of Scalp EEG"*—substantially improve the theoretical framing of the paper. The authors now clearly articulate their study's value as a head-to-head multivariate benchmark of first-order alpha amplitude against nine higher-order local cross-frequency coupling (CFC) and inter-electrode functional connectivity estimators within the same dataset and cross-validation framework. They also candidly discuss how single-trial estimator signal-to-noise ratio (SNR), continuous mouse-driven visuomotor reproduction variance, right-handedness and pooled staircase luminance, and within-analysis versus experiment-wise familywise Type-I error control shape the interpretation of their findings.

However, a rigorous audit of the revised manuscript text (`manuscript/main.tex`), the response letter (`manuscript/response_to_reviewers.tex`), the auto-generated statistical macros (`manuscript/generated_stats.tex` / `manuscript/analysis_results.json`), the newly reconstructed figures (`Figures 2–7`), and the underlying Python reconstruction scripts (`manuscript/scripts/eeg_preproc_and_feats.py` and `manuscript/scripts/run_mvpa_and_figures.py`) reveals **critical internal contradictions between the prose/macros and the newly generated figures**, as well as several **methodological discrepancies between what the Methods text states and what the analysis code actually executes**. 

Specifically, when the authors re-ran the MVPA pipeline in Python (`run_mvpa_and_figures.py`) to regenerate high-resolution Figures 3–6 and create the new summary Figure 7, the set of metrics and time windows reaching statistical significance ($p < 0.05$, $t_{\max}$ corrected) changed relative to the original MATLAB analysis. While `analysis_results.json["mvpa_summary_metrics"]` and Figures 3–7 were updated with the new Python run's outputs, **half of `analysis_results.json["mvpa"]` (`svm_results` and `svr_results`) and the corresponding Results, Discussion, Figure Legends, and Response to Reviewers text were left hardcoded to the old single-window spikes!** Consequently, the manuscript text and figure captions currently describe green significance bars and asterisks that do not exist in Figures 3–7, while failing to mention the significance bars and asterisks that *are* plotted in Figures 3–7.

Below, I first evaluate the authors' responses to my five Round 1 comments, and then detail the specific text–figure–code discrepancies that must be resolved before the manuscript can be accepted for publication.

---

## Evaluation of Responses to Round 1 Comments

### Point 1: Scalp EEG Limitations for Gamma/Connectivity & What the Study Adds Beyond Alpha Amplitude Decoding
- **Authors' Response:** Satisfactory conceptually. The authors added a dedicated subsection (*"What the Current MVPA Comparison Adds to the Literature and Signal-to-Noise Considerations"*, lines 97–100) and a new three-panel summary figure (**Figure 7**). Framing the paper as a systematic, head-to-head single-trial multivariate benchmark testing whether macroscopic alpha amplitude is merely a proxy for local alpha–gamma PAC (as posited by the OPR model; Zazio et al. 2020) or inter-areal alpha synchronization (Sadaghiani and Kleinschmidt 2016) provides a clear and defensible "take-home" message. The rewritten discussion of invasive ECoG (Voytek 2010) versus non-invasive MEG (Bonnefond and Jensen 2015) and scalp EEG myogenic/cranial attenuation constraints above $50\text{ Hz}$ (lines 107–108) is also well-reasoned.

### Point 2: Signal-to-Noise Ratio (SNR) Confounds Across EEG Metrics
- **Authors' Response:** Partially satisfactory, though purely qualitative. In the revised Discussion (line 100), the authors explicitly acknowledge that instantaneous alpha amplitude is a first-order narrowband envelope estimate of the dominant human scalp rhythm with high single-trial SNR, whereas local alpha–gamma CFC and pairwise functional connectivity are higher-order relational statistics estimated over short sliding windows ($250\text{--}375\text{ ms}$) that suffer from higher estimation variance and lower single-trial SNR.
- **Remaining Subtlety:** While this qualitative caveat is welcome and necessary, the authors did not compute any empirical proxy for estimator reliability (e.g., split-half reliability, trial-to-trial coefficient of variation, or baseline SNR across the 10 metrics) and still lean toward interpreting the null CFC/connectivity results as substantive theoretical support for BSEM over OPR in the Abstract and Conclusion. Furthermore, in `Table 1` (*Overview of predicted outcomes*), the authors state that the OPR model predicts *"At-chance classification of non-informative vs. informative cued trials by alpha amplitude alone."* This is a questionable characterization of Zazio et al. (2020): the OPR model does not deny that alpha amplitude is modulated by spatial attention cues; rather, it proposes that alpha oscillations modulate sensory responsiveness via AAC (threshold) or PAC (gain). The authors should soften `Table 1` so as not to set up a "straw-man" prediction for the OPR model.

### Point 3: Continuous $360^\circ$ Mouse-Rotation Task & Motor/Working Memory Variance
- **Authors' Response:** Satisfactory conceptually, but **undermined by a major Methods vs. Code discrepancy regarding the SVR target variable** (detailed in Critical Issue 2 below). In the Discussion (*"Task Demands, Continuous Motor Reproduction, and Visual Field Asymmetry"*, line 103), the authors thoughtfully explain how the $516.6\text{ ms}$ post-mask delay and self-paced continuous mouse rotation introduce visual working memory maintenance, visuomotor coordinate transformation, and fine motor execution variance that dilute trial-by-trial correlations ($r_z$) with pre-stimulus oscillations. However, as noted below, the manuscript text claims SVR predicted *signed* response error ($[-180^\circ, +180^\circ]$), whereas the code actually trained SVR on *absolute* response error ($|\text{error}| \in [0^\circ, 180^\circ]$).

### Point 4: Multiple Comparison Correction Across 10 Metrics and Target Sides
- **Authors' Response:** Conceptually sound in principle, **but empirically contradicted by the new Python MVPA results** (see Critical Issue 1 below). The authors explain in the Methods (line 257) and Discussion (line 100) that right-tailed $t_{\max}$ permutation testing ($10,000$ iterations, $\alpha = 0.05$) was applied across time points within each separate metric and target-side analysis, and they explicitly note that isolated transient spikes would not survive a Bonferroni threshold across the 10 metrics ($\alpha_{\text{Bonf}} = 0.005$, or $0.0025$ across 10 metrics $\times$ 2 hemifields). Ironically, in their newly executed Python pipeline (`run_mvpa_and_figures.py`), the unadjusted within-analysis $t_{\max}$ tests produced false-positive spikes in **6 out of 10 metrics for Left Targets** (`Alpha Amplitude`, `MVL_Özkurt`, `MI`, `PAC_PLV`, `FC: Correlation`, and `FC: PLV`) and **4 out of 10 metrics for Right Targets** (`Alpha Amplitude`, `MVL_Özkurt`, `MI`, and `PAC_PLV`), plus an SVR spike in `PAC_GLM`! This proliferation of 1-to-2 window spikes across half of the tested metrics in Figure 7A–B is a textbook illustration of why experiment-wise Type-I error inflation across 20 SVM and 40 SVR family tests is a serious issue—and yet the manuscript prose completely ignores these new spikes and pretends only `MVL_Özkurt` and `wPLI` had isolated spikes.

### Point 5: Visual Field Asymmetry in Behavioral Performance (Figure 2)
- **Authors' Response:** Satisfactory. In the Discussion (line 105), the authors provide two plausible, complementary explanations for why absolute orientation errors were larger on the left than the right and why spatial cueing selectively benefited left-hemifield targets: (1) all $N = 28$ participants were right-handed and operated the computer mouse with their right hand, favoring left-hemisphere (right-visual-field) sensorimotor integration during continuous visuomotor rotation; and (2) the preliminary 1-up/2-down detection staircase pooled target contrast across left and right hemifields rather than staircasing each hemifield independently, placing right-hemifield targets near a performance floor ($\sim 28^\circ\text{--}29^\circ$ error) even without an informative cue.

---

## Critical Issues Requiring Revision Before Acceptance

### 1. Severe Discrepancy Between Manuscript Text/Captions and Reconstructed Figures 3–7 (`analysis_results.json` Synchronization Bug)

In `manuscript/scripts/run_mvpa_and_figures.py` (lines 1636–1670), the script updates `analysis_results.json["mvpa_summary_metrics"]` (which drives **Figure 7**) and renders **Figures 3–6** from `results_dict`, **but only updates 5 keys inside `analysis_results.json["mvpa"]["svm_results"]` and 0 keys inside `analysis_results.json["mvpa"]["svr_results"]`**. As a result, `generated_stats.tex` and the manuscript prose (`main.tex`) are completely out of sync with Figures 3, 4, 5, 6, and 7:

1. **SVM Classification of Local CFC (Figure 3 vs. Text & Figure 7A–B):**
   - **What `main.tex` (p. 20, lines 22–23) states:** *"Among the five local alpha–gamma CFC metrics, only $\text{MVL}_{\text{Özkurt}}$ yielded above-chance SVM classification of informative versus non-informative cues, and only at isolated, transient time points: at $-875\text{ ms}$ pre-stimulus for left targets, and at $-785\text{ to }-770\text{ ms}$ and $-380\text{ ms}$ pre-stimulus for right targets (Figure 3A, fourth row; Figure 7A–B)... None of the remaining CFC measures (AAC, MI, $\text{PAC}_{\text{GLM}}$, or $\text{PAC}_{\text{PLV}}$) classified covert attention above chance at any pre-stimulus time point."*
   - **What Figure 3A and Figure 7A–B actually show:**
     - **MI** has significant green bars and inset topographies for **both Left Targets** (around $-215\text{ ms}$, $2.0\%$ of windows, blue asterisk in Fig 7A) and **Right Targets** (around $-500\text{ ms}$, $4.0\%$ of windows, red asterisk in Fig 7A).
     - **$\text{PAC}_{\text{PLV}}$** has significant green bars and inset topographies for **both Left Targets** (around $-845\text{ ms}$, $2.0\%$ of windows, blue asterisk in Fig 7A) and **Right Targets** (around $-425\text{ ms}$, $2.0\%$ of windows, red asterisk in Fig 7A).
     - **$\text{MVL}_{\text{Özkurt}}$** has **4 separate green significance bars** for Left Targets ($8.0\%$ of windows = 4/50 steps) and **5 separate green significance bars** for Right Targets ($14.0\%$ of windows = 7/50 steps)—not a single point at $-875\text{ ms}$ (which would be $1/50 = 2\%$, contradicting the $8\%$ macro `\statSVMMVLPctSigRange` right next to it!).
     - Furthermore, in **Figure 3B**, **$\text{PAC}_{\text{GLM}}$** (Left Targets) displays a green significance star on central electrode `Cz`, and **Alpha Amplitude** (Left Targets) displays a single green significance star near `CP1`/`C3`, contradicting the text's description of the significant feature weights.

2. **SVM Classification of Functional Connectivity (Figure 4 vs. Text, Caption, & Figure 7A–B):**
   - **What `main.tex` (p. 21, lines 35–45) and `response_to_reviewers.tex` state:** *"Among the four inter-electrode alpha functional connectivity metrics, only debiased weighted phase-lag index (wPLI) exceeded permuted-label chance at a single isolated time window ($-462.5\text{ ms}$ pre-stimulus for right targets; Figure 4A; Figure 7A–B)... Amplitude envelope correlation, PLV, and iCOH remained at chance throughout the entire pre-stimulus interval."* The Figure 4 caption likewise states: *"The green shaded bar indicates the single significant time window for wPLI on right-target trials ($-462.5\text{ ms}$, $p < 0.05$, $t_{\max}$ corrected)."*
   - **What Figure 4A and Figure 7A–B actually show:**
     - **wPLI** has **NO green significance bar anywhere** in Figure 4A ($0.0\%$ significant windows for both Left and Right targets in `analysis_results.json["mvpa_summary_metrics"]` and Figure 7A–B).
     - Instead, **Correlation (Left Targets)** has a green significance bar and inset connectome around $-712.5\text{ ms}$ ($8.3\%$ of windows = 1/12 steps, blue asterisk in Fig 7A), and **PLV (Left Targets)** has a green significance bar and inset connectome around $-412.5\text{ ms}$ ($8.3\%$ of windows = 1/12 steps, blue asterisk in Fig 7A)!

3. **SVR Prediction of Response Errors (Figures 5, 6, and 7C vs. Text & Captions):**
   - **What `main.tex` (pp. 22–23, lines 50–76) and `response_to_reviewers.tex` state:** SVR exceeded chance at single windows for **Alpha Amplitude** at $-640\text{ ms}$ (Non-Informative Right, with significant occipital/parietal feature weights), **AAC** at $-200\text{ ms}$ (Non-Informative Left, with significant Fz/FCz/FC1 weights), and **Correlation** at $-562.5\text{ ms}$ (Non-Informative Left, marked by a green bar in Figure 6A).
   - **What Figures 5, 6, and 7C actually show:**
     - In **Figure 5A**, neither **Alpha Amplitude** nor **AAC** has any green significance bar in any panel! Instead, the **only** green significance bar and inset topography in Figure 5A is for **$\text{PAC}_{\text{GLM}}$** under **Informative Left** (around $-800\text{ ms}$)!
     - In **Figure 6A**, there is **NO green significance bar** in any panel (contradicting both the Results text and the Figure 6 caption).
     - In **Figure 7C** (`svr_noninf_sig_indices: []`, `svr_inf_sig_indices: [4]`), the **only** asterisk is on **$\text{CFC: PAC (GLM)}$ (Informative Cue, blue bar)**—which directly contradicts the Results text (line 54: *"None of the four phase–amplitude coupling metrics ($\text{MI}$, $\text{MVL}_{\text{Özkurt}}$, $\text{PAC}_{\text{GLM}}$, or $\text{PAC}_{\text{PLV}}$) significantly predicted trial-by-trial response errors in any condition"*) and contradicts Table 1's BSEM vs. OPR evaluation!

**Action Required:** The authors must investigate why their Python re-implementation (`run_mvpa_and_figures.py`) diverged from the original MATLAB DDTBOX results (see Technical Code Issues below for the exact causes), and ensure 100% synchronization across `analysis_results.json`, `generated_stats.tex`, the Results/Discussion text, the figure captions, and Figures 3–7.

---

### 2. Methodological & Code Discrepancies Between Manuscript Text and Python Scripts

In inspecting `manuscript/scripts/eeg_preproc_and_feats.py`, `manuscript/scripts/run_mvpa_and_figures.py`, and `manuscript/scripts/run_behavior_pipeline.py` against `manuscript/main.tex`, I identified seven concrete discrepancies between what the manuscript describes and what the code actually does:

1. **SVR Target Variable: Signed Response Error (Text) vs. Absolute Response Error (Code):**
   - **In `main.tex` (p. 16, line 252):** *"For support vector regression (SVR) ($\epsilon\text{-SVR}$, cost parameter $C = 0.1$), models were trained separately within each of the four Cue Type $\times$ Target Side conditions to predict continuous single-trial **signed response errors** from the pre-stimulus spatiotemporal EEG features."*
   - **In `eeg_preproc_and_feats.py` (line 879):** `c_svr_labels = np.abs(err_deg[mask]).astype(np.float64)`.
   - Thus, the SVR models were actually trained on **absolute response error** ($|\text{Target} - \text{Response}| \in [0^\circ, 180^\circ]$), not signed response error! (Note: predicting *absolute* response error—i.e., error magnitude/precision—actually makes far more theoretical sense for testing perceptual precision under BSEM vs. OPR than predicting signed clockwise/counterclockwise bias, which has a mean of $0^\circ$. However, the Methods text must state accurately whether signed or absolute error was used.)

2. **Empirical Null Repetitions (200 in Text vs. 100 in Code) and Null Variance Inflation (`tmax_classifier_accuracies.m`):**
   - **In `main.tex` (p. 16, line 254) and `generated_stats.tex` (`\statCVRepeatsNull{20}`, `\statCVTotalNull{200}`):** States that empirical chance distributions were estimated using 10-fold CV repeated **20 times (200 null evaluations)** per participant.
   - **In `run_mvpa_and_figures.py` (lines 310, 405–421, 1402):** Runs only **10 repeats $\times$ 10 folds = 100 null evaluations** (`perm_acc` shape `(n_steps, 10, 10)`).
   - **Crucially:** In `run_mvpa_and_figures.py` (line 1402), the script averages across all 100 null evaluations (`all_perm_mean = np.mean(all_perm_reps, axis=(2, 3))`) before computing `diff_scores = all_real - all_perm_mean`. In DDTBOX (`tmax_classifier_accuracies.m`, as the comment on line 1403 even notes!), group testing either compares against a single drawn permutation repetition (10-fold average) or accounts for null variance. Averaging 100 null shuffles collapses the variance of `all_perm_mean` around $50.0\%$ (or $0.0$ for SVR) by a factor of $\sqrt{10}$, artificially shrinking the standard error of the difference (`stder_obs`) and inflating $t$-statistics across all metrics! This explains why the gray error bars on the permutation chance traces in Figures 3–6 are nearly invisible and why so many spurious single-window $t_{\max}$ spikes appeared in Figure 3, Figure 4, Figure 5, and Figure 7.

3. **Haufe Feature Weight Normalization Bug (`/ 20.0` vs. `/ 100.0`):**
   - In `run_mvpa_and_figures.py` (lines 426 and 584), `sum_abs_fw` accumulates `abs(ap)` across `10 repeats * 10 folds = 100` cross-validation iterations, plus across `win_len` time samples via `np.mean(...)`, and then divides by `20.0` (`ch_abs[ch] = np.mean(sum_abs_fw[ch * win_len : (ch + 1) * win_len]) / 20.0`). While the subsequent within-step $z$-scoring (`fw_z[s] = (ch_abs - m_ch) / sd_ch`) makes a constant scalar divisor cancel out for single steps, this shows the code was adapted from an earlier 2-repeat $\times$ 10-fold draft. More importantly, taking `abs(ap)` on every single fold *before* averaging across folds (`sum_abs_fw[p] += abs(ap)`) turns zero-mean noise weights into a positive Rician/folded-normal bias proportional to feature variance, which differs from averaging signed Haufe patterns across folds before taking the magnitude.

4. **Figure 3A & 5A Inset Topography Permutations ($2,000$ in Code vs. $10,000$ in Text) and Max-2-Insets Truncation:**
   - **In `main.tex` (p. 17, line 259):** States that feature weights for significant windows were tested using $10,000$ permutations.
   - **In `run_mvpa_and_figures.py` (lines 998–1001):** Uses `n_perm=2000` (`mxt_perm1_phipson_smyth(avg_run_z, n_perm=2000, ...)`), and `if run_idx < 2:` silently drops inset topographies for the 3rd, 4th, and 5th significant windows (e.g., in Alpha Amplitude Right Targets and $\text{MVL}_{\text{Özkurt}}$ Left and Right Targets in Figure 3A, where green bars appear without inset heads).

5. **Behavioral Analysis Trial Inclusion (All Completed Trials vs. EEG-Clean Trials):**
   - **In `main.tex` (p. 11, line 184):** *"All epochs were then visually inspected to reject trials (from both EEG and behavioral analyses) contaminated by amplifier blocking, residual blinks, muscle bursts, or skin potentials."*
   - **In `run_behavior_pipeline.py` (lines 67–88):** Loads raw behavioral `.mat` files (`{sub}_Orient_Exo.mat`), removes only the 20 practice trials (and post-343 truncation for subject `009`), and computes mean `abs_err` and `resp_err` across **all behavioral trials without applying the EEG/HEOG/visual-inspection trial rejections** (`keep1`, `keep2`, `keep3`, `keep4` in `eeg_preproc_and_feats.py`).

6. **Preprocessing Order & Baseline Correction Omission in `eeg_preproc_and_feats.py`:**
   - **In `main.tex` (p. 11, line 182):** States that epochs were baseline-corrected using the $200\text{ ms}$ interval following fixation onset, rejected at $\pm 1000\,\mu\text{V}$, ocular-corrected via Gratton et al. (1983), **followed by a second $200\text{ ms}$ post-fixation baseline subtraction**, and rejected at $\pm 500\,\mu\text{V}$.
   - **In `eeg_preproc_and_feats.py` (lines 768–826):** Applies a 4th-order Butterworth IIR bandpass filter (`sosfiltfilt`, rather than EEGLAB's FIR filter stated in line 182) and performs **neither the first nor the second $200\text{ ms}$ post-fixation baseline subtraction** before checking the $\pm 1000\,\mu\text{V}$ and $\pm 500\,\mu\text{V}$ thresholds! Because EEG epochs can have DC/slow-drift offsets prior to baseline subtraction, omitting baseline subtraction alters which trials exceed $\pm 500\,\mu\text{V}$ or $\pm 1000\,\mu\text{V}$.

7. **Figure 2 Panel Structure (3 Panels A, B, C in Image vs. 2 Panels A, B in Text & Caption):**
   - **In `figure2_behavior.png` (generated by `run_behavior_pipeline.py`, lines 175–232):** The figure has **three panels**: **A. Absolute Orientation Error**, **B. Directional Bias**, and **C. Staircase Threshold** (a histogram of participant target luminance thresholds, Mean = 54.0).
   - **In `main.tex` (Figure 2 caption, lines 8–10, and Results text, lines 267–268) and `response_to_reviewers.tex` (line 101):** Only **Panel A** and **Panel B** are described; **Panel C** is never mentioned in the Figure 2 caption or anywhere in the Results text! Furthermore, the Figure 2A caption states: *"($*p < 0.05$, FDR-corrected $t_{\max}$ permutation tests; $\text{ns} = \text{not significant}$)"*, but `figure2_behavior.png` does not actually draw any `*` or `ns` significance brackets on Panel A!

---

## Minor Presentation & Editorial Issues

1. **Color Mismatch Between Figures 3–6 and Their Captions:**
   - The captions for Figures 3, 4, 5, and 6 state: *"Black traces denote empirical accuracy [or $r_z$]; **blue traces** denote permuted-label chance accuracy [or null $r_z$]."*
   - In the newly rendered Figures 3–6 (`run_mvpa_and_figures.py`, lines 1015 and 1181), the permuted-label chance traces are plotted in **gray** (`color="#7F7F7F"`), not blue, and there is no legend drawn in Panel A. Please update either the caption ("gray dashed traces") or the trace color, and add a small legend to Panel A.
2. **Vertical Dashed Line at $0\text{ ms}$ in Figure 3–6 Captions:**
   - The Figure 3 caption states: *"The vertical dashed line marks target onset ($0\text{ ms}$)."* Because `xlim` ends at `0.0` and the right spine is hidden (`ax.spines["right"].set_visible(False)`), no vertical dashed line is visible at $0\text{ ms}$ in Figures 3–6.
3. **Broken Citation Author Name (Typo in `main.tex`):**
   - On page 16 (`main.tex`, line 250), `(Haufe et al. 2014; Grootswagers et al. 2017;itchie et al. 2019)` is missing the `R` in `Ritchie et al. 2019`.
4. **FOOOF Frequency Range Discrepancy (`main.tex` vs. `eeg_preproc_and_feats.py`):**
   - `main.tex` (line 190) states FOOOF was fit across `3.9 to 48.8 Hz` (`\statFOOOFFreqRange`), whereas `eeg_preproc_and_feats.py` (line 597) fits across `2.0 <= freqs[k] <= 50.0` Hz (`1.95 Hz` resolution from `nfft=512` at $1000\text{ Hz}$).

---

## Summary of Required Revisions

To make `CerCor-2022-00354R1` acceptable for publication in *Cerebral Cortex*, the authors should:
1. **Fix the MVPA Null Comparison & Synchronize All Text, Macros, and Figures:** Correct the permutation null baseline in `run_mvpa_and_figures.py` so that it matches DDTBOX (`20 repeats` and proper null variance handling rather than averaging across 100 null shuffles), update all fields of `analysis_results.json["mvpa"]` and `analysis_results.json["mvpa_summary_metrics"]` from the same analysis run, and ensure that every significant window and asterisk shown in **Figures 3–7** matches the **Results text**, **Discussion**, **Figure Captions**, and **Response to Reviewers**.
2. **Reconcile Methods vs. Code Discrepancies:** Explicitly clarify whether SVR predicted *absolute* or *signed* response error (aligning `main.tex` line 252 and `eeg_preproc_and_feats.py` line 879), ensure behavioral trial rejection matches the Methods description (or state accurately how behavioral trials were filtered), and ensure baseline subtraction in `eeg_preproc_and_feats.py` matches the EEG Preprocessing text.
3. **Harmonize Figure 2 and Figures 3–6 Captions:** Add `*` / `ns` annotations to Figure 2A and describe **Figure 2C** (staircase luminance distribution) in the Figure 2 caption and Results text; update the Figures 3–6 captions to match the gray dashed permutation traces; and fix the typo in `Ritchie et al. 2019`.
