# Old World calibration points for lexical dating (clepsydra)

Companion to `calibrations_oldworld.tsv`, which has 153 rows: 75 node calibrations, 8 stem/originate calibrations and 70 tip (attested-language) ages. It covers Europe, Africa, the Near East, the Caucasus, South/East/Southeast Asia and Siberia. Compiled 2026-10-10.

## 1. Conventions

- **All ages are years BP with 1950 CE = 0.** Most sources use a different "present", so every value was shifted, and the column `original_spec_as_published` keeps the value as it was published:
  - Holman et al. 2011, Rama & Wichmann 2020, IE-CoR/Heggarty 2023, Chang et al. 2015, Kolipakam et al. 2018, Savelyev & Robbeets 2020, Kitchen et al. 2009 and Bouckaert et al. 2012 use present ≈ 2000, so I subtracted 50. For IE-CoR I checked this against two tips: OHG at 1170 BP gives 830 CE (the Tatian translation), and Gothic at 1630 BP gives 370 CE. Bouckaert is "≈2000" by inference from its Latin = 2050.
  - Lee & Hasegawa 2011 use present ≈ 2010, inferred from "Tokugawa 1603 CE = 407 YBP". I subtracted 60.
  - Zhang et al. 2019 state present = 1950, so no shift.
  - Sagart et al. 2019, Lee 2015 and Maurits et al. 2020 don't state a base year. I left them unshifted or approximated them, and the notes say which.
- **min / max.** `age_min_BP` is the youngest age the node can have: the split happened *no later than* this. `age_max_BP` is the oldest it can have: the split happened *no earlier than* this. Watch the terminology in Maurits et al. 2020: their "upper bound for the divergence" means my **max** age.
- **Attestation logic** follows Maurits et al. 2020:
  - A dated text in language L gives a **min** for every clade that contains L's lineage, together with a sister.
  - A text in a *unified* proto-language, or a loan into it, gives a **max** on that proto-language's break-up.
  - A migration into new territory gives a **max** on the split it causes.
  - Point dates derived from events are usually one of these bounds mis-specified as a point. `bound_type` flags this ("point (max-leaning)").
- **Grades.**
  - A: a strong historical or epigraphic anchor, used with the right logic.
  - B: archaeological or loanword inference with an explicit argument, or a historical event loosely tied to the split.
  - C: contested, a correlation between an archaeological culture and a language, a political event used as a proxy, or a prior far narrower than its evidence.
  - U: the value was verified but I found no justification for it (Zhang 2019 Karenic, Ngwi, Pumi and Qiang–rGyalrong; Kolipakam's Brahui).
- **`suggested_distribution`** gives either the published prior, shifted to 1950, or a better-posed alternative. Where only a single date exists I added the rule of thumb from Maurits et al. 2020: the 95% width is max(age/2, 1500 yr), with Normal sd = width/3.92.
- **Glottocodes** were looked up by name in glottolog-cldf `languages.csv` (master branch, fetched 2026-10-10). Tip glottocodes come from IE-CoR. "Dardic" isn't a Glottolog clade, so the code given is for Eastern Dardic.

## 2. Sources actually read (and how)

| Source | What I extracted | Access |
|---|---|---|
| Holman et al. 2011, *Curr. Anthropol.* 52:841, doi:10.1086/662127 | All 52 calibrations, with dates, type (A/H/E) and justifications. I kept the 40 in scope and skipped the American, Oceanic and Australian ones and Inuit. | Full PDF (MPG PuRe) |
| Rama & Wichmann 2020, *PLoS ONE* 15:e0236522, doi:10.1371/journal.pone.0236522 | Table 1: 30 non-overlapping groups, which are a subset of Holman's. Table 7: leave-one-out results for 8 ancient IE tips. | PMC full text |
| Heggarty et al. 2023, *Science* 381:eabg0818, doi:10.1126/science.abg0818 | The priors for all 52 ancient tips (Normal or offset-lognormal). | IE-CoR CLDF `languages.csv` (lexibank/iecor). I didn't read the SM text. |
| Bouckaert et al. 2012, *Science* 337:957, doi:10.1126/science.1219669 | 20 tip priors and 13 node priors, read directly from the published BEAST XML. | phlorest/bouckaert_et_al2012 |
| Chang et al. 2015, *Language* 91:194, doi:10.1353/lan.2015.0005 | 13 uniform tip priors and 8 ancestry constraints. There are no node-age priors. | phlorest/chang_et_al2015 `ieo.xml` (narrow dataset) |
| Gray & Atkinson 2003, *Nature* 426:435, doi:10.1038/nature02029 | 14 constraints (11 nodes plus the Hittite and Tocharian tips) are mentioned, but the SI table was not accessible. | Accepted MS (Oxford ORA). See §5. |
| Kitchen et al. 2009, *Proc R Soc B* 276:2703, doi:10.1098/rspb.2009.0408 | 5 tip dates, 4 stem windows and the root window. | PMC text and archived XML. The two disagree, see §4. |
| Kolipakam et al. 2018, *R Soc Open Sci* 5:171504, doi:10.1098/rsos.171504 | Table 2 (8 priors), cross-checked against the XML. | PMC and phlorest XML |
| Sagart et al. 2019, *PNAS* 116:10317, doi:10.1073/pnas.1817972116 | 4 tips from the text, plus a Sinitic node prior that appears only in the XML. | PMC and phlorest XML |
| Zhang et al. 2019, *Nature* 569:112, doi:10.1038/s41586-019-1153-z | 8 priors from the XML. I didn't read the justifications in the SI. | phlorest XML |
| Lee & Hasegawa 2011, doi:10.1098/rspb.2011.0518; 2013, doi:10.1371/journal.pone.0062243; Lee 2015, doi:10.1371/journal.pone.0128448 | Japonic, Ainu and Koreanic priors. | PMC (and XML for Ainu) |
| Honkola et al. 2013, doi:10.1111/jeb.12107, and Maurits et al. 2020, *JoLE* 5:17, doi:10.1093/jole/lzz009 | The original and revised Uralic calibrations, plus the methodological critique. | Maurits PDF (OUP). I know Honkola's values only as Maurits reports them. |
| Savelyev & Robbeets 2020, *JoLE* 5:39, doi:10.1093/jole/lzz010 | Old Turkic and Cuman tips. | PDF (OUP) |
| Grollemund et al. 2015, doi:10.1073/pnas.1503793112; Koile et al. 2022, doi:10.1073/pnas.2112853119 | The 4 Bantu calibrations and Koile's re-specification. | PMC |
| Serva & Petroni 2008, *EPL* 81:68005, doi:10.1209/0295-5075/81/68005 | Two point calibrations: Italian–French 1600 y and Icelandic–Norwegian 1100 y. | arXiv |

I checked all DOIs above against Crossref.

## 3. Quality: what is solid and what isn't

**Strongest material (grade A).**

- **Ancient-language tips** with securely dated texts. These include Hittite, Mycenaean, Old Persian, Latin, Ancient Greek, Gothic, OCS, OHG, Old English, Old Norse, Classical Armenian, Old Irish, Old Chinese, Old Japanese, Middle Korean, Old Turkic, Akkadian, Ugaritic, Old Khmer, Old Georgian and Old Javanese.
- **Hard minimums from attestation.** Examples:
  - Tamil-Brahmi (254 BCE) gives a min for South Dravidian I.
  - Akkadian gives a min for Semitic.
  - Old Iranian texts give a min for Iranian.
  - Mitanni Indo-Aryan gives a min for Indo-Iranian. I proposed this one myself, so check its exact date before use.
- **Hard maximums from settlement of new land.** Examples are Iceland (874 CE), the Cape (1652), Malta's Arab conquest (870 CE) and the Madagascar migration.

**Methodologically sound inferential calibrations (B).**

- Maurits et al. 2020 derive Finnic, Saami and Permic from loanword stratigraphy and argue each explicitly, from at least two independent lines of evidence. These are the best-documented node calibrations in the set.
- Aikio's Saami date, used by Holman, matches Maurits' Saami calibration almost exactly: 1700 BP both ways.

**Weak (C). Use with caution, and don't use these as held-out "truth".**

- **Archaeological cultures equated with proto-languages.** This covers Parpola's equations for Indo-Iranian, Indo-Aryan, Iranian and Dardic, the Benue-Congo and Bantu calibrations from Shum Laka, Obobogo and Urewe, Southern Nilotic, and the IE minus Anatolian/Tocharian node. Maurits et al. 2020 single out Grollemund 2015 as the kind of practice to avoid. Holman's Indo-Aryan (3900) and Dardic (3550) calibrations had two of the worst ASJP discrepancies of all 52, at −95% and −90%. Only Czech–Slovak, at −119%, was worse.
- **Political events used as proxies.** Examples: 476 CE for Italo-Western Romance, the Turk Khaganate for Common Turkic, and the Kipchak empire.
- **Implausibly narrow priors in Bouckaert 2012.** Lithuanian–Latvian has sd 25, Brythonic sd 25, and West Germanic sd 25 on a node that probably doesn't match the event. Northwest Germanic is N(1875, 67), which Maurits criticises explicitly. Balto-Slavic has two priors stacked on it.
- **Loanword maximums mis-specified as points.** Hmong-Mien is dated by the loan of Chinese "money" (5th c. BCE), which logically gives only a max.

## 4. Conflicts: same node, different values (all in 1950 BP)

| Node or tip | Values |
|---|---|
| Romance | Holman 1679 (Dacia 271 CE, point). Serva & Petroni: Italian–French = 1600 "time distance". Bouckaert calibrates Latin+Romance instead, N(1950, 135). |
| Brythonic | Holman 1400 (point, from a 5th–7th c. range). Bouckaert N(1500, 25). |
| Northwest/West Germanic | Bouckaert NW Germanic N(1825, 67) and West Germanic N(1500, 25). Holman's "Germanic" is 2050 but is really NW Germanic, since it has no Gothic. The three can't all be right: West Germanic can't be 1500 if NW Germanic is 1825 and English–Frisian (Holman) is 1500. |
| Slavic | Holman 1400 (point). Bouckaert offset-lognormal with min 1150 and mean ≈1450. |
| Indo-Aryan, Iranian, Indo-Iranian | Holman points of 3850, 3850 and 4350 (Parpola). Bouckaert gives only minimums: 2100, 2550 and 2950. |
| Saami | Holman 1700 and Maurits N(1700, 153). These agree. |
| Permic | Honkola U(1100, 1300), shifted. Maurits reversed-lognormal with 95% in 661–1263. |
| Samoyedic, Finno-Saamic | Honkola U(2000, 2200) and 2000–3000. Maurits dropped both. |
| Ethiosemitic | Holman 2400 (calibration). Kitchen estimated 2850 (posterior). |
| Sinitic | Holman 1950 (Qin–Han). Sagart XML U(2000, 2200), base year unclear. Zhang calibrates a different node, Mandarin + Old Chinese, at N(2700, 150). |
| Old Chinese tip | Sagart U(2300, 2800). Zhang N(2500, 100). These are compatible. |
| Ainu, Sakhalin | Paper text says N(500, ~51). The published XML says N(460, 51). |
| Semitic tips and stems | Kitchen's paper uses uniform 1000-year windows, e.g. Hebrew origin 3200–4200. The phlorest XML has Normal(·, 250) priors on *leaf heights* with different means (Hebrew 3800). I followed the paper. Kitchen's "Biblical Aramaic = 1800 YBP" also looks too young. |
| Gothic tip | IE-CoR N(1580, 10), Bouckaert N(1600, 25), Chang U(1575, 1625). |
| Hittite tip | IE-CoR N(3400, 125), Bouckaert N(3400, 125), Chang U(3250, 3450). |
| Vedic tip | IE-CoR offset-lognormal with offset 2950. Bouckaert N(2950, 100). Chang U(2950, 3450). |
| Middle Korean tip | Lee 2015 N(500, 100) truncated, with no base year stated. |
| Bantu | Grollemund fixes East Bantu at a point of 2500 (archaeology). Koile 2022 drops that calibration and turns the others into lognormals. |

The alternatives for every IE-CoR tip are in the `alt_values_other_studies` column.

## 5. Not verified or not found (flagged rather than guessed)

- **Gray & Atkinson 2003 SI table.** I couldn't access the list of the 14 constraints. Bouckaert 2012 largely inherits them, and those values were verified from its XML.
- **Supplementary justifications I didn't read:** Heggarty 2023 SM (justification for each tip), Bouckaert 2012 SI, Zhang 2019 SI, Kolipakam SI §1 (Brahui), and Koile 2022's lognormal parameters (SI Fig. S6).
- **Kassian et al. 2021** (*Linguistics* 59:949, doi:10.1515/ling-2020-0060). The publisher blocked access, so none of its calibrations are included.
- **Starostin-school dates.** These are glottochronological *estimates*, not calibrations, so I didn't include them.
- **Families with no calibrations found:**
  - Tai-Kadai, Nakh-Daghestanian and Berber: no published Bayesian dating with calibrations, and no suitable Holman points.
  - Kartvelian: none either, apart from the Old Georgian tip.
  - Austroasiatic: none found beyond the Old Khmer tip. Old Mon and Middle Vietnamese (the 1651 de Rhodes dictionary) are possible tips, but I didn't verify them.
  - Mongolic and Tungusic: only the shallow Holman points.
- **Proposed by me, not taken from a dating study:** Iceland–Norway max, Afrikaans–Dutch max, Mitanni min for Indo-Iranian, Arabic-dialects max, and the Old Khmer, Old Javanese and Old Georgian tips. Their sources are general historical or epigraphic references. Old Javanese's earliest date was reportedly revised (Kesongo 743 CE), and I saw that only in a search summary.

## 6. Recommendations for clepsydra

1. **Held-out benchmarks (use c).**
   - Prefer grade A tips together with grade A/B *bounded* nodes: Finnic, Saami, Permic, South Dravidian I min, Semitic min, Maltese–Maghrebi, Ma'anyan–Malagasy, East Slavic, the Scandinavian split, Romani dispersal, Oromo and Swahili.
   - Score a benchmark as "inside the interval", not as error against a point.
   - Follow Rama & Wichmann and use non-overlapping groups, so that nested calibrations aren't double-counted.
2. **Cross-family rate priors (use b).** Draw from as many families as possible. The table is heavily weighted towards Indo-European (81 of 153 rows, 52 of them IE-CoR tips). If a rate prior pools these naively, the result will effectively be an IE rate.
3. **Encode the bound direction explicitly.** Every Holman point marked "max-leaning" should become a one-sided or asymmetric prior. Maurits et al. 2020 show that this kind of re-specification moved the Proto-Uralic estimate by about 850 years and widened its HPD by about 3000 years.
4. **Datum.** All Holman-derived points are midpoints of source ranges and use present = 2000. The ranges themselves are in Holman's "Source" text and in this file's `event` column. Use them instead of the midpoints where possible.
