# Cerebral Cortex — Round 2 Editorial Consultation & Conference Call Transcript

**Manuscript ID:** CerCor-2022-00354R1  
**Title:** *"Connecting Covert Attention and Visual Perception to the Spatiotemporal Dynamics of Alpha Band Activity, Cross-Frequency Coupling (CFC), and Functional Connectivity using Multivariate Pattern Analysis (MVPA)"*  
**Authors:** Sarah S. Sheldon, Alona Fyshe, and Kyle E. Mathewson (University of Alberta)  
**Date of Editorial Conference:** November 9, 2026  
**Format:** Zoom Video/Audio Editorial Consultation (Recorded & Transcribed for Editorial File)  

**Participants:**
- **Dr. Pasko Rakic, MD, PhD (Editor-in-Chief, *Cerebral Cortex*)** — Chairing the consultation; focused on scientific rigor, broad neurobiological significance, journal standards, and actionable guidance for the authors.
- **Reviewer 1 (Visual Attention, Cortical Oscillations & Psychophysics Specialist)** — Focused on theoretical models of alpha oscillations (BSEM vs. OPR), psychophysical task dynamics ($360^\circ$ continuous orientation reproduction vs. 2-AFC), estimator signal-to-noise ratio (SNR), and experiment-wise familywise error control.
- **Reviewer 2 (M/EEG Signal Processing, Cross-Frequency Coupling & MVPA Methods Specialist)** — Focused on electrophysiological accuracy (ECoG vs. MEG vs. scalp EEG), time-frequency/spectral parameterization (FOOOF, sideband-preserving PAC filters, debiased wPLI), DDTBOX cross-validation mechanics, and figure-to-text fidelity.

---

## Segment 1: Opening Remarks & Assessment of Round 1 $\rightarrow$ Round 2 Improvements

**Dr. Pasko Rakic (Editor-in-Chief):**  
Good afternoon to both of you, and thank you for making the time for this live editorial consultation on manuscript `CerCor-2022-00354R1` from Kyle Mathewson's group at the University of Alberta—Sarah Sheldon, Alona Fyshe, and Kyle Mathewson. I have read both of your exceptionally thorough Round 2 reports, and I have also spent the morning inspecting their revised manuscript (`manuscript/main.tex`), their point-by-point response letter (`manuscript/response_to_reviewers.tex`), the newly rendered figures, and their open-source reproducibility repository (`OrientExoML`).  
*[Inner thought: It is rare to get two reviewers who not only read every line of the prose and response letter, but actually pull down the authors' GitHub/OSF code and cross-check the JSON macros against the rendered PNGs. We have a fascinating situation here: conceptually and visually, the authors did everything we asked for in Round 1—and yet their new automated Python pipeline created a glaring split-brain mismatch between the figures and the text. I want to make sure we acknowledge the authors' genuine hard work first before we dissect the technical synchronization bugs.]*  
Before we dive into the statistical and code-synchronization anomalies that you both uncovered, let's establish where we stand on the conceptual, structural, and pedagogical revisions we requested in Round 1. Reviewer 1, in the first round you had five substantive concerns—about what a scalp EEG study can truly add given gamma and connectivity limitations, single-trial SNR confounds, motor variance in the $360^\circ$ mouse-rotation task, multiple comparison correction across 10 metrics, and the left-versus-right visual field asymmetry in Figure 2. How satisfied are you with how the authors addressed the *conceptual* side of those five points?

**Reviewer 1:**  
Thank you, Pasko. On the conceptual and theoretical front, I have to give the authors a lot of credit. They took every single one of my Round 1 critiques seriously and wrote three substantial new subsections in the Discussion—*"What the Current MVPA Comparison Adds to the Literature and Signal-to-Noise Considerations"*, *"Task Demands, Continuous Motor Reproduction, and Visual Field Asymmetry"*, and *"Implications for the OPR Model and Electrophysiological Scope of Scalp EEG"*.  
*[Inner thought: Honestly, in Round 1 I was leaning toward rejection or a very skeptical major revision because I felt the paper was dressed-up alpha lateralization with nine noisy null results tacked on. But their new framing—combined with Reviewer 2's idea for Figure 7—actually convinced me of the paper's value as a head-to-head methodological and theoretical benchmark.]*  
Specifically, framing the paper as the first systematic, head-to-head single-trial multivariate benchmark comparing first-order alpha amplitude against five local CFC estimators and four volume-conduction-mitigated functional connectivity estimators within the exact same participants ($N = 28$), task epochs, and cross-validation framework gives the paper a clear, defensible "take-home" message. Their discussion of single-trial estimator SNR (contrasting narrowband alpha envelopes against $250\text{--}375\text{ ms}$ sliding-window relational statistics), their explanation of how the $516.6\text{ ms}$ post-mask blank delay and self-paced mouse rotation introduce working memory and visuomotor noise, and their two-part explanation of the right-visual-field behavioral advantage (right-handed mouse control favoring left-hemisphere visuomotor integration, plus a pooled detection staircase that put the right visual field near a $\sim 28^\circ$ floor) are all thoughtful and convincing.

**Reviewer 2:**  
I completely agree with Reviewer 1 on the conceptual and stylistic overhaul. In Round 1, I had a long list of complaints about colloquial phrasing ("the idea is", "As can be seen in..."), the tangled structure of the Introduction, missing methodological parameters, a muddled paragraph in the Discussion that conflated subdural ECoG theta–high-gamma PAC (Voytek 2010) with MEG alpha–high-gamma PAC (Bonnefond and Jensen 2015) and "surface EEG", and the fact that Figures 3 through 6 were practically illegible postage stamps.  
*[Inner thought: They really went above and beyond on the writing and figure layout. Most authors just bump the font size by 2 points in Illustrator and call it a day; Sheldon and Mathewson actually built an entire end-to-end programmatic pipeline with a 3-panel summary figure (Figure 7). That's exactly why it's such a shame that their Python re-implementation introduced four signal-processing/statistical bugs and left `analysis_results.json` half-updated!]*  
They fixed every single stylistic and pedagogical point I raised. The Abstract is crisp; the Introduction cleanly separates the background literature on BSEM (Samaha et al. 2020), OPR (Zazio et al. 2020), and network connectivity (Sadaghiani and Kleinschmidt 2016) from their own hypotheses in the final paragraph; the Methods section now gives the exact mathematical definitions of signed versus absolute orientation error, the $25\,\mu\text{V}$ split-half sliding-window HEOG step threshold for overt saccades, the FOOOF $1/f$ versus Gaussian peak parameterization, the grand-average IAF statistics ($M = 10.6\text{ Hz}$, $\text{SD} = 0.5\text{ Hz}$, range $9.6\text{--}11.4\text{ Hz}$), and the sideband-preservation rationale ($\pm 0.4 \times f_c$, giving an effective $15\text{--}50\text{ Hz}$ passband) for the high-frequency gamma filter (Aru et al. 2015; Seymour et al. 2017) along with an honest caveat about upper-beta ($15\text{--}30\text{ Hz}$) leakage. And visually, the new 300 DPI Figures 3–6 and the new **Figure 7** summary comparison are night-and-day better than Round 1.

**Dr. Pasko Rakic (Editor-in-Chief):**  
And they also formatted all references strictly to *Cerebral Cortex* style and embedded all figures and tables directly in the double-spaced manuscript as requested.  
*[Inner thought: Good—both reviewers recognize that the core science and scholarly presentation are now at a high level. That means we are NOT looking at a rejection; we are looking at a targeted revision to fix the pipeline/text synchronization and code discrepancies.]*  
Now, let us turn to the heart of why we cannot accept the manuscript in its current state (`CerCor-2022-00354R1`): the internal contradictions between the manuscript text/macros and the newly generated Figures 3–7.

---

## Segment 2: The "Split-Brain" Discrepancy Between Manuscript Text (`generated_stats.tex`) and Reconstructed Figures 3–7

**Dr. Pasko Rakic (Editor-in-Chief):**  
Both of you independently caught something extraordinary. When the authors built their open-science Python pipeline (`manuscript/scripts/run_mvpa_and_figures.py`) to regenerate high-resolution Figures 3–6 and construct the new Figure 7, the script updated `analysis_results.json["mvpa_summary_metrics"]` (which drives Figure 7) and rendered Figures 3–6 from the new Python run, **but left almost all of `analysis_results.json["mvpa"]["svm_results"]` and `["svr_results"]` untouched from their earlier MATLAB run**. Reviewer 2, walk us through how that happened in the code and what it did to the paper.

**Reviewer 2:**  
Yes, if you look at `manuscript/scripts/run_mvpa_and_figures.py` at lines 1636–1670, the script calls `api.handle_request({"action": "update_stats", ...})`. Inside that dictionary, it updates `mvpa_summary_metrics` completely, and inside `mvpa["svm_results"]` it updates only five keys: `alpha_pct_sig_range` (`"90\\%--96\\%"`), `alpha_pct_sig_min` (`90`), `alpha_left_sig_windows`, `alpha_right_sig_windows`, and `mvl_pct_sig_range` (`"8\\%--14\\%"`). It **never updates** `mvl_left_sig_ms`, `mvl_right_sig_window_ms`, `mvl_right_sig_point_ms`, `wpli_pct_sig_range`, `wpli_right_sig_ms`, or **any** of the keys inside `mvpa["svr_results"]` (`alpha_noninf_right_sig_ms`, `aac_noninf_left_sig_ms`, `corr_noninf_left_sig_ms`)!  
*[Inner thought: As someone who writes analysis pipelines, I know exactly how this happens: you write a JSON-to-LaTeX macro generator (`generate_latex_stats.py`) around the original manuscript's results structure, and then when you re-run the Python port to generate the figures, you only wire up the summary arrays and forget that the prose has hardcoded sentences wrapped around the old single-window macros.]*  
Because `generate_latex_stats.py` reads `analysis_results.json` to create `generated_stats.tex`, the compiled manuscript (`main.tex`) and the `response_to_reviewers.tex` letter are literally describing a different set of statistical results than what is plotted right in front of the reader's eyes in Figures 3, 4, 5, 6, and 7!

**Reviewer 1:**  
Let's list the three specific places where the prose and captions clash head-on with Figures 3–7, because it is jarring for any careful reader:
1. **First, Local CFC SVM Classification (Figure 3 vs. Text & Figure 7A–B):**  
   In `main.tex` (lines 289–290), the text states that among the five local CFC metrics, *"only $\text{MVL}_{\text{\"Ozkurt}}$ yielded above-chance SVM classification... at $-875\text{ ms}$ pre-stimulus for left targets, and at $-785\text{ to }-770\text{ ms}$ and $-380\text{ ms}$ pre-stimulus for right targets"* and that *"None of the remaining CFC measures (AAC, MI, $\text{PAC}_{\text{GLM}}$, or $\text{PAC}_{\text{PLV}}$) classified covert attention above chance at any pre-stimulus time point."*  
   Yet when you look at **Figure 3A** and **Figure 7A–B**:
   - **MI** has green significance bars and inset topographies for **both Left Targets** (at $-215\text{ ms}$, $2.0\%$ of windows, blue asterisk in Fig 7A) and **Right Targets** (at $-500\text{ to }-485\text{ ms}$, $4.0\%$ of windows, red asterisk in Fig 7A)!
   - **$\text{PAC}_{\text{PLV}}$** has green significance bars and inset topographies for **both Left Targets** (at $-845\text{ ms}$, $2.0\%$ of windows, blue asterisk in Fig 7A) and **Right Targets** (at $-425\text{ ms}$, $2.0\%$ of windows, red asterisk in Fig 7A)!
   - **$\text{MVL}_{\text{\"Ozkurt}}$** has **4 separate significant windows** for Left Targets ($-875, -695, -560, \text{ and }-395\text{ ms}$, which equals $4/50 = 8.0\%$ of windows) and **5 separate clusters (7 steps)** for Right Targets ($-800, -725, -470, -440, \text{ and }-380\text{ to }-350\text{ ms}$, which equals $7/50 = 14.0\%$ of windows). Even worse, in the Figure 7B caption (line 352), the macro `\statSVMMVLPctSigRange` expanded to `$8\%--14\%$`, while the prose describes it as a "single-window" effect at $-875\text{ ms}$ (which would be $2\%$)!
   - And in **Figure 3B**, `PAC_GLM` (Left Targets) has a green significance star on central electrode `Cz` that is never mentioned in the text, while `Alpha Amplitude` (Left Targets) only shows one green star near `CP1`/`C3` instead of the broad centroparietal/parieto-occipital cluster described in the text.

2. **Second, Functional Connectivity SVM Classification (Figure 4 vs. Text, Caption, & Figure 7A–B):**  
   In `main.tex` (lines 302–312) and the Figure 4 caption, the authors state that **only wPLI** exceeded chance at a single window ($-462.5\text{ ms}$ on right-target trials) and that *"The green shaded bar indicates the single significant time window for wPLI on right-target trials ($-462.5\text{ ms}$)"*, while Correlation, PLV, and iCOH remained at chance.  
   Yet in **Figure 4A** and **Figure 7A–B**, **wPLI has NO green bar at all ($0\%$ on both Left and Right)**! Instead, **Correlation (Left Targets, $-712.5\text{ ms}$, $8.3\%$)** and **PLV (Left Targets, $-412.5\text{ ms}$, $8.3\%$)** have green significance bars, 3D connectome insets, and blue asterisks in Figure 7A–B!

3. **Third, SVR Prediction of Response Errors (Figures 5, 6, & 7C vs. Text & Captions):**  
   In `main.tex` (lines 317–343), the text and captions state that SVR exceeded chance for **Alpha Amplitude** at $-640\text{ ms}$ (Non-Informative Right), **AAC** at $-200\text{ ms}$ (Non-Informative Left), and **Correlation** at $-562.5\text{ ms}$ (Non-Informative Left, explicitly referenced in the Figure 6 caption), and states unequivocally (line 321): *"None of the four phase–amplitude coupling metrics ($\text{MI}$, $\text{MVL}_{\text{\"Ozkurt}}$, $\text{PAC}_{\text{GLM}}$, or $\text{PAC}_{\text{PLV}}$) significantly predicted trial-by-trial response errors in any condition."*  
   Yet in **Figure 5A**, **Figure 6A**, and **Figure 7C**, Alpha Amplitude, AAC, and Correlation have **zero** green bars and **zero** asterisks! Instead, the **only** green significance bar and asterisk in the entire SVR analysis is on **$\text{PAC}_{\text{GLM}}$ under Informative Left** (at $-800\text{ ms}$, blue bar with asterisk in Figure 7C)!  
   *[Inner thought: And notice the theoretical irony: in Round 1 I asked them about familywise Type-I error inflation across 10 metrics. Their new Python run—because of the null variance deflation bug we're about to discuss—sprayed 1-to-2 window false-positive spikes across 6 of the 10 metrics in SVM and onto PAC_GLM in SVR, which visually undermines their own text claiming only MVL_Özkurt and wPLI had isolated spikes!]*

**Dr. Pasko Rakic (Editor-in-Chief):**  
That is a crystal-clear diagnosis. Obviously, we cannot publish a paper in *Cerebral Cortex* where the figure captions tell the reader to look at a green bar on wPLI in Figure 4A and Correlation in Figure 6A that isn't there, while ignoring green bars on MI, $\text{PAC}_{\text{PLV}}$, Correlation, PLV, and $\text{PAC}_{\text{GLM}}$ that *are* there. Now, why did the Python pipeline (`eeg_preproc_and_feats.py` and `run_mvpa_and_figures.py`) produce different significance masks than the original MATLAB pipeline (`repos/MVPA_CuedOrientTask_paper/Analysis_Code/`)? Let's walk through the methodological and code discrepancies you both discovered.

---

## Segment 3: Methodological & Code Discrepancies Uncovered in the Repository Audit

**Dr. Pasko Rakic (Editor-in-Chief):**  
Reviewer 2, in your report you pinpointed the exact mathematical and signal-processing differences between the authors' original MATLAB codebase and their new Python scripts, and Reviewer 1, you independently found the same null-variance bug, the SVR signed-versus-absolute error discrepancy, the Haufe weight accumulation bug, and the behavioral trial-filtering mismatch. Let's go through these systematically so we can give the authors an exact, actionable engineering checklist.

### 3.1. The Group-Level Permutation Baseline Variance Deflation Bug (`run_mvpa_and_figures.py`)

**Reviewer 2:**  
The single biggest reason why Figures 3A, 4A, 5A, and 7A–B sprouted so many new isolated green significance bars in the Python run is in `manuscript/scripts/run_mvpa_and_figures.py` at lines 310, 421, 1402–1404, and 1444–1445.  
First, the manuscript text (line 254) and `generated_stats.tex` (`\statCVRepeatsNull{20}`, `\statCVTotalNull{200}`) state that the empirical chance distribution was estimated using 10-fold CV repeated **20 times (200 null evaluations)** per participant, which matches `permut_rep = 20` in the original MATLAB scripts (`Run_ConnDecoding_v1.m` and `Run_SVR_v2.m`). In `run_mvpa_and_figures.py`, they only ran **10 repeats $\times$ 10 folds = 100 null evaluations**.  
Second—and much more importantly—look at how the group-level paired $t_{\max}$ test compares empirical accuracy against permuted-label null accuracy. In DDTBOX (`Analyse_Decoding_Results_v2.m` line 223: `analysis_settings.perm_disp.drawmode = 2`, which calls `tmax_classifier_accuracies.m` lines 32–34), for each participant, DDTBOX draws **one random cross-validation repetition** (an average of 10 folds) from that participant's permuted-label distribution (`all_subj_perm_acc_reps_draw`) to compare against the real cross-validated accuracy, so that the empirical score and the chance baseline have comparable estimation variance across participants.  
In contrast, in `run_mvpa_and_figures.py` (line 1402: `all_perm_mean = np.mean(all_perm_reps, axis=(2, 3))`), the Python script **averages across all 100 null evaluations** per participant before computing `diff_scores = all_real - all_perm_mean`!  
*[Inner thought: Even the comment on line 1403 of their own script says `# Draw random repetition per subject as in tmax_classifier_accuracies.m (default: mean or drawn perm)`, which shows the programmer knew what DDTBOX did, but took the shortcut of subtracting `all_perm_mean`!]*  
When you average 100 null evaluations per subject, the variance of `all_perm_mean` around $50\%$ (or $0.0$ for SVR) shrinks by a factor of $10$ (standard deviation shrinks by $\sqrt{10}$). That artificially deflates the standard error of `diff_scores` (`stder_obs`) in the denominator of the paired $t$-statistic, inflating $t$-values across all time steps and creating false-positive single-window spikes in `MI`, `PAC_PLV`, `Correlation`, `PLV`, and `PAC_GLM`. It also explains why the gray error bars (`se_perm`) on the permuted chance traces in Figures 3–6 are almost invisibly flat compared to the empirical error bars!

**Reviewer 1:**  
Spot on. When the authors fix `run_mvpa_and_figures.py` to run 20 null repeats (200 evaluations) and use DDTBOX's `drawmode = 2` (or match the exact variance structure of `tmax_classifier_accuracies.m`), two good things will happen at once: the permutation SEM error bars in Figures 3–6 will look realistic again, and the spurious single-window false positives that contradicted their prose will largely disappear.

### 3.2. SVR Target Variable: Signed Response Error (Text) vs. Absolute Response Error (Code)

**Reviewer 1:**  
The second major discrepancy is in `manuscript/main.tex` on page 16 (line 252):  
*"For support vector regression (SVR) ($\epsilon\text{-SVR}$, cost parameter $C = 0.1$), models were trained separately within each of the four Cue Type $\times$ Target Side conditions to predict continuous single-trial **signed response errors** from the pre-stimulus spatiotemporal EEG features."*  
When I checked `manuscript/scripts/eeg_preproc_and_feats.py` (line 879), the code actually computes:  
`c_svr_labels = np.abs(err_deg[mask]).astype(np.float64)`  
And as Reviewer 2 verified in the original MATLAB repository (`repos/MVPA_CuedOrientTask_paper/Analysis_Code/PrepData_DDTBOX.m`, lines 115–118), the original MATLAB code *also* used `abs(errordeg)`!  
*[Inner thought: Thank goodness the code used `abs(errordeg)`! If they had actually trained a linear SVR on pre-target cue-locked EEG before the random target orientation was even chosen to predict clockwise ($+$) versus counterclockwise ($-$) rotational bias, the entire SVR analysis would have been scientifically nonsensical.]*

**Reviewer 2:**  
Exactly! That was a huge relief when I checked `PrepData_DDTBOX.m`. Because the target orientation ($0^\circ\text{--}360^\circ$) is chosen randomly on every trial *after* the pre-stimulus interval, pre-stimulus oscillatory activity before target onset cannot possibly know whether a participant will err clockwise or counterclockwise—and more importantly, the theoretical question differentiating BSEM (Samaha et al. 2020) from OPR (Zazio et al. 2020) is whether pre-stimulus oscillations predict **objective perceptual precision** (i.e., error magnitude, $|\text{Target} - \text{Response}|$). So the code in both MATLAB and Python did the right scientific analysis (`abs(err_deg)`), and the authors simply wrote "signed response errors" by mistake in line 252 of `main.tex`. They just need to change "signed response errors" to **"absolute response errors (error magnitude)"** in the Methods, Results, and Figure 5–7 captions.

### 3.3. Biased Standard wPLI vs. Vinck's Debiased Squared wPLI (`dwPLI`)

**Reviewer 2:**  
Next, let's look at why the wPLI functional connectivity trace changed between MATLAB and Python. In `main.tex` (Table 2, line 216; Methods, line 247; Results, line 302), the authors explicitly state that they computed the **debiased weighted phase-lag index (wPLI)** of Vinck et al. (2011). Indeed, in their original MATLAB repository (`conn_wPLI.m` line 64, called by `PrepConnData_DDTBOX.m` line 21), they computed Vinck's debiased squared wPLI estimator:
$$\text{dwPLI} = \frac{\left(\sum_{t} \text{Im}(X_t)\right)^2 - \sum_{t} \text{Im}(X_t)^2}{\left(\sum_{t} |\text{Im}(X_t)|\right)^2 - \sum_{t} \text{Im}(X_t)^2}$$
However, in `manuscript/scripts/eeg_preproc_and_feats.py` (line 551), the Python port wrote:
```python
wpli_out[ij, pair_idx, tr] = abs(wpli_num) / wpli_den
```
which is the **biased standard wPLI** ($|\sum \text{Im}(X_t)| / \sum |\text{Im}(X_t)|$), omitting the sample-size debiasing terms $-\sum \text{Im}(X_t)^2$! Replacing those five lines in `eeg_preproc_and_feats.py` with the exact Vinck et al. (2011) debiased formula from `conn_wPLI.m` is trivial and will restore exact mathematical equivalence between the Python pipeline and the manuscript Methods.

### 3.4. Missing Pre- and Post-EMCP $200\text{ ms}$ Baseline Subtractions and Filter Type in `eeg_preproc_and_feats.py`

**Reviewer 1:**  
We also both caught a discrepancy in the EEG preprocessing script (`manuscript/scripts/eeg_preproc_and_feats.py`, lines 768–826) versus the EEG Preprocessing text (`main.tex`, lines 182–184):
- `main.tex` (line 182) states that continuous EEG was filtered with EEGLAB's FIR filter, segmented, **baseline-corrected using the $200\text{ ms}$ interval following fixation onset**, thresholded at $\pm 1000\,\mu\text{V}$, ocular-corrected via Gratton et al. (1983) EMCP, **followed by a second $200\text{ ms}$ post-fixation baseline subtraction**, and thresholded at $\pm 500\,\mu\text{V}$.
- In `eeg_preproc_and_feats.py` (lines 778–825), both the pre-EMCP and post-EMCP $200\text{ ms}$ baseline subtractions (`pop_rmbase` in `Preproc_OrientWheel_Exo_ML.m` lines 92 and 171) are omitted before checking the $\pm 1000\,\mu\text{V}$ and $\pm 500\,\mu\text{V}$ voltage rejection thresholds! Because un-baselined EEG epochs have DC offsets, omitting those two baseline subtractions changes which trials exceed $\pm 500\,\mu\text{V}$ or $\pm 1000\,\mu\text{V}$.

### 3.5. Behavioral Pipeline Trial Filtering (`run_behavior_pipeline.py` vs. `main.tex` line 184)

**Reviewer 1:**  
Similarly, in `main.tex` (line 184), the authors state: *"All epochs were then visually inspected to reject trials (from both EEG and behavioral analyses) contaminated by amplifier blocking, residual blinks, muscle bursts, or skin potentials."*  
Yet in `manuscript/scripts/run_behavior_pipeline.py` (lines 67–88), the script computes mean absolute error and signed error across all completed behavioral trials in `{sub}_Orient_Exo.mat` (after removing the 20 practice trials and subject `009`'s post-343 trials) without applying the EEG artifact rejection masks (`keep1`–`keep4`). The authors should either apply the exact same trial-rejection mask in `run_behavior_pipeline.py` or clarify in `main.tex` (line 184) that behavioral ANOVAs were conducted on all completed behavioral trials while MVPA used the artifact-free EEG trials (matching `perm_rbANOVA_OrientExo.m` in `CuedOrientTask_paper`).

### 3.6. Haufe Feature Weight Normalization (`/ 20.0` vs. `/ 100.0`) & Inset Topography Permutation Count (`2000` vs. `10000`)

**Reviewer 1:**  
Two more quick code-versus-text details in `run_mvpa_and_figures.py`:
1. In `ddtbox_svm_subject` (line 426) and `ddtbox_svr_subject` (line 584), `sum_abs_fw` accumulates `abs(ap)` over `10 repeats * 10 folds = 100` folds and then divides by `20.0` instead of `100.0` (`ch_abs[ch] = np.mean(sum_abs_fw[...]) / 20.0`), and takes `abs(ap)` per fold before averaging rather than averaging the Haufe activation pattern across folds before taking the magnitude or testing against permuted feature weights.
2. In `main.tex` (line 259), the text states that feature weights for significant windows were tested using $10,000$ permutations (`\statMVPAPerms`), whereas `run_mvpa_and_figures.py` (line 999) calls `mxt_perm1_phipson_smyth(avg_run_z, n_perm=2000, ...)` and caps the number of plotted inset topographies at two per panel (`if run_idx < 2:` on line 1001). Changing `n_perm=2000` to `n_perm=10000` takes less than a second in Numba and makes the code match the text 100%.

### 3.7. FOOOF Frequency Range (`3.9 to 48.8 Hz` in Text vs. `2.0 to 50.0 Hz` in Code)

**Reviewer 2:**  
And finally, on FOOOF: `main.tex` (line 190) and `analysis_results.json` report that FOOOF was fit across `3.9 to 48.8 Hz` (`\statFOOOFFreqRange`) with peak width limits `[1, 5] Hz`, max peaks `7`, min peak height `0.22`, and peak threshold `0.20`. Meanwhile, in `eeg_preproc_and_feats.py` (line 597), the custom Numba routine fits across `2.0 <= freqs[k] <= 50.0` Hz (which at $\Delta f = 1000/512 = 1.953\text{ Hz}$ corresponds to bins from $1.95\text{ to }48.83\text{ Hz}$, or $3.91\text{ to }48.83\text{ Hz}$ if excluding the lowest bin). The authors should ensure the frequency range and settings in `eeg_preproc_and_feats.py` and `main.tex` are described consistently.

---

## Segment 4: Theoretical & Presentation Refinements (Table 1, Figure 2, Figures 3–6 Captions, and Typos)

**Dr. Pasko Rakic (Editor-in-Chief):**  
*[Inner thought: That list of seven code/method items is remarkably clear and constructive. Because the authors have already built a modular pipeline (`OrientExoML`), fixing those items and re-running `generate_latex_stats.py` will bring the entire paper into lockstep.]*  
Let's now review the remaining theoretical and presentation refinements that you both noted in `Table 1`, `Figure 2`, and `Figures 3–6`.

**Reviewer 1:**  
First, in **Table 1** (*Overview of predicted outcomes that support the BSEM or the OPR model*, page 6, lines 108–136), under the right-hand column (**Oscillation-Based Probability of Response (OPR) Model**), Prediction 1 currently states:
> *"At-chance classification of non-informative vs. informative cued trials by alpha amplitude alone AND Above-chance classification of non-informative vs. informative cued trials by phase–amplitude coupling (PAC) OR amplitude–amplitude coupling (AAC)."*

As I noted in my review, claiming that the OPR model (Zazio et al. 2020) predicts *at-chance* classification of spatial attention cues by alpha amplitude alone is overly strict—almost a straw-man prediction. Zazio et al. (2020) do not deny that spatial attention modulates macroscopic alpha amplitude; rather, they propose that alpha amplitude/phase influences visual perception *via* local AAC (modulating threshold) or PAC (modulating response gain). So if OPR were the primary mechanism linking attention to perception, AAC or PAC should *also* strongly classify attentional states (comparable to or mediating alpha amplitude) and, critically, **Prediction 2** is the true decisive test: under OPR, PAC (or AAC) must predict trial-by-trial perceptual response errors, whereas under BSEM, pre-stimulus alpha amplitude shifts baseline excitability (classifying attentional cues) *without* predicting objective orientation discrimination precision on non-informative trials. Softening Prediction 1 in the OPR column of Table 1 so that it doesn't claim OPR requires zero alpha-amplitude modulation during spatial cueing will make their theoretical framing airtight.

**Reviewer 2:**  
I agree completely with Reviewer 1 on Table 1. Now let's look at the visual presentation issues in **Figure 2** and **Figures 3–6**:
1. **Figure 2 (`figure2_behavior.png` vs. `main.tex` lines 268–278):**  
   - **Unreferenced Panel C:** When the authors regenerated `figure2_behavior.png` using `run_behavior_pipeline.py` (lines 175–232), they included **three panels**: **A. Absolute Orientation Error**, **B. Directional Bias**, and **C. Staircase Threshold** (a histogram of participant target grayscale values, showing `Mean = 54.0`). However, the Results text (lines 268–269), the Figure 2 caption (lines 276–277), and the Response to Reviewers (line 101) only describe **Panel A** and **Panel B**—**Panel C is never mentioned anywhere!** Honestly, Panel C is quite a nice addition (especially since Reviewer 1 and the authors discuss the staircase threshold in the Discussion!), so they should simply add a sentence describing **Figure 2C** to the Figure 2 caption and reference `(Figure 2C)` in the Methods/Results where the staircase grayscale distribution ($M = 54.0, \text{SD} = 22.0$) is reported.
   - **Missing `*` and `ns` Significance Annotations on Panel A:** The Figure 2 caption explicitly states: *"Informative spatial cues significantly reduced error magnitude for left-hemifield targets ($*p < 0.05$, FDR-corrected $t_{\max}$ permutation tests; $\text{ns} = \text{not significant}$)."* Yet in `figure2_behavior.png`, `run_behavior_pipeline.py` did not draw the `*` bracket over the Left Target bars or the `ns` bracket over the Right Target bars! Adding those two significance brackets in `run_behavior_pipeline.py` will make Figure 2A match its caption.

2. **Figures 3–6 Captions and Visual Polish (`main.tex` lines 298, 311, 329, 342):**  
   - **"Blue traces" vs. Gray Dashed Traces for Permuted Chance:** All four captions (Figures 3, 4, 5, and 6) still say *"Black traces denote empirical accuracy [or $r_z$]; **blue traces** denote permuted-label chance accuracy [or null $r_z$]."* In the newly rendered Figures 3–6, the permuted-label chance traces are plotted in **gray dashed lines (`#7F7F7F`)**, while in Figures 5 and 6 **blue (`#0072B2`)** is used for the empirical Informative Cue trace and **magenta/red (`#D55E00` / `#CC79A7`)** is used for the Non-Informative Cue trace! Calling the null trace "blue" in Figures 5 and 6 is actively confusing. The captions for Figures 3–6 must be updated to accurately describe the trace colors (e.g., gray dashed traces for permuted-label chance, and specifying the exact colors for Informative vs. Non-Informative in Figures 5–6), and a compact legend should be added to Panel A of Figures 3–6.
   - **"Vertical dashed line marks target onset (0 ms)":** The captions for Figures 3–6 mention a vertical dashed line at $0\text{ ms}$, which is not visible because `xlim` ends at `0` and the right axis spine is hidden. Either draw a subtle vertical dashed line at $x = 0\text{ ms}$ or remove that phrase from the captions (and mention the horizontal dashed reference line at $50\%$ chance or $r_z = 0.0$ instead).
   - **Inset Topography Nose Clipping & Max-2-Insets Cap:** In `run_mvpa_and_figures.py` (line 1018), lowering the inset y-coordinate slightly (e.g., from `0.56` to `0.50` or `0.52`) will prevent the tip of the nose on the circular topoplots in Figures 3A, 4A, and 5A from being clipped by the top subplot border.

3. **Typo in `main.tex` Line 250:**  
   - Fix the broken citation `(Haufe et al. 2014; Grootswagers et al. 2017;itchie et al. 2019)` $\rightarrow$ `Ritchie et al. 2019`.

---

## Segment 5: Consensus Editorial Decision & Action Plan

**Dr. Pasko Rakic (Editor-in-Chief):**  
Let's synthesize our consensus decision. Both of you agree that the scientific premise, the experimental dataset ($N = 28$), the multivariate benchmarking design across 10 oscillatory/CFC/connectivity metrics, the expanded theoretical Discussion, and the visual redesign (including the new Figure 7 summary) are now strong and well-suited for *Cerebral Cortex*.  
*[Inner thought: If we send a crystal-clear, itemized decision letter that gives the authors the exact line numbers in `main.tex`, `analysis_results.json`, `eeg_preproc_and_feats.py`, `run_mvpa_and_figures.py`, and `run_behavior_pipeline.py`, they can execute these fixes cleanly and rapidly.]*  
However, we are all in complete agreement that the manuscript **cannot be accepted until the text–figure–code synchronization errors and the four underlying Python translation bugs are fully resolved**. Specifically:
1. Every number, significant time window, percentage of significant windows, and asterisk in **Figures 3, 4, 5, 6, and 7** must be **100% identical** to what is stored in `analysis_results.json` / `generated_stats.tex` and what is described in the **Results**, **Discussion**, **Figure Captions**, and **Response to Reviewers**.
2. The Python scripts (`eeg_preproc_and_feats.py`, `run_mvpa_and_figures.py`, and `run_behavior_pipeline.py`) must be corrected so that they faithfully implement the methods stated in the manuscript (specifically: DDTBOX `drawmode = 2` / matched-variance permutation baseline with 20 repeats rather than averaging 100 null shuffles; Vinck et al. 2011 debiased squared wPLI; $200\text{ ms}$ baseline subtractions before and after Gratton EMCP; `n_perm = 10000` for inset feature weight tests; and harmonized FOOOF frequency range and behavioral trial filtering).
3. The Methods text (line 252) must be corrected from **"signed response errors"** to **"absolute response errors"** to match what both the MATLAB and Python SVR pipelines actually predicted.
4. **Table 1**, **Figure 2** (adding `*`/`ns` brackets on Panel A and describing Panel C in the caption/text), **Figures 3–6** (fixing caption trace color descriptions, adding a Panel A legend, and preventing inset nose clipping), and the **`Ritchie et al. 2019`** typo must be cleaned up.

Reviewer 1, Reviewer 2—do you both endorse a decision of **Minor-to-Moderate Revision (Conditional Acceptance pending Technical & Text–Figure Synchronization Verification)** with this exact checklist?

**Reviewer 1:**  
Yes, I fully endorse that decision. Once the null-variance bug is fixed and the text, JSON macros, and Figures 3–7 are 100% synchronized, this will be a very solid contribution to the oscillations and visual attention literature.

**Reviewer 2:**  
I completely concur. Giving the authors the exact file names and line numbers will make it straightforward for them to fix the Python port, regenerate the figures and LaTeX macros from a single unified run, and return a clean, publication-ready manuscript.

**Dr. Pasko Rakic (Editor-in-Chief):**  
Wonderful. Thank you both for an exemplary, collegial, and deeply rigorous review process. I will draft the official Round 2 Decision Letter incorporating our unified checklist and both of your complete review reports right away.

---
*End of Editorial Conference Call Transcript — Recorded for Manuscript File CerCor-2022-00354R1.*
