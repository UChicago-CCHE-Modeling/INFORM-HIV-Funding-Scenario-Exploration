Subject: Revised analysis using your respondent-level r summaries

Hi Anna,

Thank you for the summary file (pre_mtg_survey_summary_gtz_20261019). I replaced the numbers built from multiplying medians with your respondent-level effective reductions: the median of r_i is now the central point (PrEP and ART channels), and the surrogate is queried directly at P1 (ART only), P2 (PrEP only) and P3 (both). I used only the median, minimum and maximum from the file. The review PDF is attached.

Two caveats are stated in the report. First, because only the marginal summaries of r_i are available, the pairing of PrEP and ART values by respondent is unknown, so the lower and upper points combine the marginal minima and maxima; they are drawn as bars, not a joint range, and no single respondent need have reported them. Second, the survey wording is open to interpretation, so the PDF also shows r = q as a sensitivity. q is described as a gross disruption, not an observed funding reduction, and the funding-space panels are hypothetical what-ifs. ART adherence remains a provisional proxy for viral suppression.

**No individual-level respondent data have been opened, committed to Git, or processed through Claude, Codex, or any other AI system.**

If you and AFC approve, I can next finalize a local Python or R workflow that reads the respondent CSV directly on my computer and generates one short PDF per respondent. I would first test the workflow using synthetic records. The real file and all derived outputs would remain in a private local folder: they would not be committed to Git, uploaded to an AI system, or sent anywhere by the script. AFC could then relay each PDF to the corresponding respondent and ask whether their stated reduction already accounted for access by other means. I have not run this workflow on the real respondent file.

For this review, I would especially appreciate confirmation that the revised interpretation is what you intended and whether you would like me to proceed with the local respondent-level validation. It would also be helpful to confirm the source and intended interpretation of the illustrative Medicaid-loss value used in the supplementary policy panels.

Best,
Raff
