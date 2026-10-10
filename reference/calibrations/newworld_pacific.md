# Calibration points: Americas, Pacific, Australia–New Guinea, creoles

Companion to `calibrations_newworld_pacific.tsv` (214 rows, 20 columns, tab-separated, `.` = missing).
Compiled 2026-10-10 for clepsydra. Every date is either taken from a cited source or, where this
compilation did arithmetic (unit conversion, glottochronological formula on published percentages),
the row says so. Unverifiable items are graded **U**.

## 1. How to read the table

| column | meaning |
|---|---|
| `id` | prefix = block: `ARW` Arawan, `ASJP` ASJP benchmark ages, `HOL` Holman et al. 2011 calibration points, `AN` Austronesian/Pacific, `NA` North America, `MA` Mesoamerica, `SA` South America, `AU` Australia, `PG` New Guinea, `CR` creoles and colonial splits |
| `type` | `node_split`, `tip_attestation`, `root`, `settlement_bound` |
| `age_BP_min/max/point` | years BP, **1950 = 0** (CE year *y* → 1950 − *y*; post-1950 events are negative) |
| `bound` | `min` = node **no younger than** (split happened **no later than** this); `max` = node **no older than** (split **no earlier than**); `range`; `tip` |
| `suggested_distribution` | clepsydra types from `docs/calibrations.md` (`uniform_node`, `normal_node`, `offset_lognormal_node`, `soft_minimum`, `soft_maximum`, `hard_*`, `fixed_tip`) |
| `grade` | **A** strong historical record; **B** archaeological / well-argued inference; **C** contested, or a published *estimate* usable only as benchmark; **U** unverified |

Counts: A = 53, B = 41, C = 107, U = 13. By type: 71 node_split, 63 root, 58 tip_attestation, 22 settlement_bound.

### Time-basis traps (important)
* **Holman et al. 2011** (ASJP and HOL rows) give dates as years before **2000**; converted here (−50).
* **Gray et al. 2009 / Greenhill et al. 2010** calibration ranges (AN rows citing them) are *as published*, on an implicit ~2000–2010 present (Favorlang/Siraya c.1650–1662 entered as 346–384 BP). Shift by ~50–60 y if you need strict 1950 basis; negligible against their widths except for tips.
* **Greenhill et al. 2023 Uto-Aztecan** give Molina 1555 as 445 BP (2000 basis); the TSV uses 395. Their node offsets may share the 2000 basis (~50 y; negligible).

### Root rows are benchmarks, not calibrations
All 63 `root` rows are **published estimates** (ASJP, glottochronology, Bayesian posteriors). They are for
held-out benchmarking (use c) or plausibility windows. Using them as calibrations in a lexical-dating
analysis is circular. The ASJP block (39 rows) is the only uniform-method source covering nearly all
South American families, including Arawan.

## 2. Recommended core set (directly usable, grade A/B)

* **Tips (A):** Old Tupi/Tupinambá (Anchieta 1595), Old Guaraní (Ruiz de Montoya 1639/40), Colonial Quechua (Santo Tomás 1560; Third Lima Council 1584; Huarochirí c.1608), Colonial Aymara (Bertonio 1603/1612), Puquina (Oré 1607), Mochica (Carrera 1644), Muisca (Lugo 1619), Mapudungun (Valdivia 1606), Classical Nahuatl (Olmos 1547; Molina 1555/1571), Colonial Yucatec, Kaqchikel (Annals 1571–1604), Colonial Zapotec (Córdova 1578), Purépecha (Gilberti 1558–59), Massachusett (Eliot 1663), Laurentian (Cartier 1535–36), Huron (Sagard 1632), Taíno, Old Island Carib (Breton 1665), Favorlang, Siraya, Old Malay (682 CE), Old Javanese, Malagasy (1603), Tagalog (1593), Sydney language (Dawes 1790–91), Arawá (Chandless 1867). Tip ages matter only if the lexical data are actually coded from those sources (Ferraz Gerardi et al. 2023 *rejected* Tupinambá/Old Guaraní as tips because each conflates sources spanning >100 years).
* **Nodes (A/B):** Proto-Polynesian ≤ Tonga settlement 2838 BP (U/Th, Burley et al. 2012); East Polynesian short chronology (Wilmshurst et al. 2011); Māori vs Tahitic ≥ NZ settlement c.1280 CE; Proto-Oceanic / Lapita onset 3470–3250 cal BP (Denham et al. 2012); Uto-Aztecan subgroup offsets (Greenhill et al. 2023); Wati (Pama-Nyungan) gamma(offset 3000) (Bouckaert et al. 2018); Michael et al. 2023's 12 archaeological Arawakan calibrations (pre-Taíno Uniform 2445–2800 and pre-Island Carib Uniform 1490–1690 are true node calibrations; the others are sampled-ancestor constraints split here into crown-max/stem-min); Inuit/Thule ≥1200 CE; Cholan split by 400 CE (epigraphy); Cham, Tsat–N. Roglai (historical); Pitkern/Norfuk 1856; Sranan/Saramaccan 1651–c.1712; Mauritian/Seychellois 1721–1770.

## 3. Regional notes (condensed)

### Austronesian / Pacific (AN, 36 rows; plus HOL)
* Gray, Drummond & Greenhill 2009 used 14 calibrations. The Science SOM (Table S3) was not retrievable (403); the ranges were recovered from a secondary copy quoting Table S3 (Kaiping's `burstclock` repository) and partly cross-checked against Greenhill et al. 2010 (5 "robust" calibrations) and Greenhill et al. 2017 (Normal priors). Treat as verified-secondary.
* **Conflicts:** (i) East Polynesia — Gray's 1150–1800 BP (long chronology) vs Wilmshurst 2011 short chronology implying Proto-EP c.700–930 BP; Gray's range biases Polynesian ages upward. Holman's East Polynesian 1050 (pre-2011) is likewise superseded. (ii) Hawaii — Wilmshurst AD 1219–1266 vs Athens et al. 2014 AD 940–1130, non-overlapping. (iii) Proto-Oceanic 3200–3600 (2009/2010) vs 3000–3300 (Gray et al. 2011). (iv) PMP 3600–4500 vs N. Luzon initial pottery 5430–4290 cal BP (Cochrane et al. 2021); Gray 2010 dropped PMP as contestable. (v) Chamic 1800–2500 criticised by Maurits et al. 2020 (attests a state, not a proto-language). (vi) Marianas c.1500 BCE (Carson) vs a later Rieth & Athens model (abstract unread).
* Settlement dates are **one-sided**: settlement of island X is a *maximum* for X-internal diversification, and a *minimum* for the split of X's language from its source only if contact ceased. Gray often entered them as two-sided ranges.
* Greenhill 2017's PAn N(5200,300) is Gray 2009's posterior → benchmark only.

### North America (NA, 21 rows) and Mesoamerica (MA, 16 rows)
* Only Uto-Aztecan has a published, parameterised Bayesian calibration set (Greenhill et al. 2023, *Language* 99(1), doi:10.1353/lan.2023.0006): NUA, Takic, Tübatulabal, Numic, Nahuan offset-lognormals + Classical Nahuatl tip; root 4095 y (95% HPD 3258–5025), rejecting the ~8900 y "ancient origin". Figures read from the OSF preprint; check against the published supplement.
* No dated Bayesian trees found for Mayan, Otomanguean, Mixe-Zoque, Algonquian, Iroquoian, Siouan, Athabaskan, Eskimo-Aleut — their roots are glottochronology/expert inference (grade C/U).
* Thule: classical ~1000 CE replaced by ≥1200 CE spread (Friesen & Arnold 2008). Apachean: only firm anchor is the Querechos of 1541. Mayan: Kaufman c.2200 BCE vs ASJP 2220 y — a 2-ky disagreement. Mixe-Zoque/Olmec link and Epi-Olmec decipherment contested.

### South America (SA, 54 rows; plus ARW, ASJP)
* **Arawakan:** Walker & Ribeiro 2011 is *undated* (no calibrations). The usable set is **Michael et al. 2023, *Interface Focus* 13:20220049** (12 archaeological calibrations, ~100 historical/ethnonym tips with adjusted BP and glottocodes in the ESM — worth importing wholesale). Possible inconsistency: Santa Paula cited at 1530±30 BP in text but calibrated Normal(900,60).
* **Tupí-Guaraní:** Ferraz Gerardi et al. 2023 (PLOS ONE e0272226) used **no internal calibrations** — a uniform root prior (bounds unstated) plus tip-age adjustment for old collections; results TG crown 1762 BP (847–2740), Awetí split 2603, Mawé split 3312. Conflict: Tupiguarani pottery in the Paraná basin 2090–1740 cal BP and Rodrigues' glottochronological Proto-TG c.2500 BP. Michael et al. 2015 is undated.
* **Quechua:** Inca radiocarbon starts (Marsh et al. 2017) predate the ethnohistoric 1463/1471; using them as maxima for Ecuadorian Kichwa or Santiago del Estero Quichua presumes Inca introduction (contested).
* **Gaps:** no published absolute dates (Bayesian or glottochronological) located for Panoan, Tacanan, Cariban, Jê/Macro-Jê, Chocoan, Barbacoan, Yanomami, Mataguayan, Uru-Chipaya beyond ASJP. Tukanoan figures (2500–3000 BP) exist only in conference slides (U).
* Ceramic–language equations (Saladoid/Barrancoid ↔ Arawakan, Amazonian Polychrome ↔ Tupí or Arawakan, Taquara-Itararé ↔ Southern Jê) are the weak link; Michael et al. 2023 is the only work stating explicit criteria. Use soft priors. Taylor & Rouse 1955 (IJAL, doi:10.1086/464319) is the early glottochronology-vs-Saladoid attempt for the West Indies.

### Australia and New Guinea (AU 7, PG 2)
* Bouckaert, Bowern & Atkinson 2018 used **one** node calibration (Wati: gamma offset 3000, α=2, β=359; breakup ≤1.5 ka) plus tip dates; root 5671 (4455–6966) scales ~linearly with it → never use AU-02 (root) together with AU-01 (Wati). Holman's Pama-Nyungan 4500 (Evans & Jones 1997) is an alternative, contested anchor.
* Tasmania: Bass Strait closure ≥11,960 (prob. ~12,400) y, but Bowern 2012 finds Tasmanian languages not demonstrably related to the mainland, so the bound rarely attaches to a node.
* TNG: no dated phylogeny; Pawley & Hammarström 2018 "perhaps 10,000–6,000 years ago" via Kuk agriculture (phases themselves disputed).

### Creoles and colonial splits (CR, 17 rows)
Creoles are not tree-like (lexifier + substrates + continued contact). Use only as short-depth benchmarks or rate checks, never as calibrations inside non-creole trees. Colonial dialect splits (Brazilian/European Portuguese, Québec French, American English, Rioplatense, Afrikaans/Dutch) use founding dates only as **maximum** bounds.

## 4. Arawan and its neighbours: what is and is not known about time depth

### Composition (Glottolog 5, `araw1282`)
Arawá† (`arua1263`); Madi (`jama1261`: Jarawara `jara1276`, Banawá `bana1307`, Eastern Jamamadí `east2800`); Madihá (`madi1263`: Kulina/Madija `culi1244`, Deni `deni1241` incl. Western Jamamadí `paui1240`); Paumarí (`paum1247`); Suruwahá (`suru1263`). Madi + Madihá = Madi-Madihá (`madi1262`). This follows Dienst 2008, *LIAMES* 8:61–67 (doi:10.20396/liames.v8i1.1471), who found evidence only for Madi-Madihá as a subgroup, with Arawá probably close to it; Paumarí and Suruwahá unplaced. Jolkesky 2016 gives a different internal split (Suruwahá vs Madi-Deni-Paumarí).

### What exists
1. **No Bayesian, no archaeological, no historical calibration for any Arawan node has been published.** No Arawan node has an associated datable event: no colonial-era documentation before the 19th century, no epigraphy, no accepted archaeological correlate.
2. **ASJP (Holman et al. 2011, Table 7):** Arauan, 17 pairs, 4 of 5 subgroups, mean similarity 28.88 → **1764 y before 2000 ≈ 1714 BP** (ARW-003). ASJP errors average 29% of age, so ±~500 y at 1 sd.
3. **Dixon's cognate percentages** (Dixon 2004, *Anthropological Linguistics* 46:1–83, p.3; reproduced by Dienst 2008), 115-item list:

   |  | Arawá | Madi | Madihá | Paumarí |
   |---|---|---|---|---|
   | Madi | 79 | | | |
   | Madihá | 72 | 61 | | |
   | Paumarí | 54 | 50 | 43 | |
   | Suruwahá | 53 | 37 | 34 | 34 |

   Arawá's figures rest on Chandless's 52-item list and are less reliable. **Not a published date**, but applying the Swadesh–Lees formula to the deepest pairs (34%) gives ~2.5 ky (r = 0.805/ky) to ~3.6 ky (r = 0.86/ky) (ARW-004, grade U; arithmetic done here). The list is not a Swadesh list.
4. **Jolkesky 2016** (UnB PhD thesis, §5.3.6) gives a scenario: "proto-proto-Arawa" and proto-Harakmbut-Katukina lived in contiguous areas on the right bank of the Solimões until the middle of the second half of the 2nd millennium BCE. Proto-Arawa itself then emerged within a Central Amazonian interaction sphere (with Arawakan, Tupí, Mura, Yanomami, Maku, Jivaroan, Taruma and Nadahup contacts). Jolkesky explicitly rejects absolute dating of proto-languages. ARW-005 turns his scenario into a weak soft maximum of ~3200 BP (C).
5. **KAHD** (Ferraz Gerardi, Aragon & Reichert 2022, *JOHD* doi:10.5334/johd.80): the family "appears to be quite shallow" from low cognate diversity. No date given. Its UPGMA tree reproduces Dienst's classification.
6. **Tip ages:** Arawá known only from Chandless 1867 (83 BP; ARW-001, grade A). Ehrenreich's Purus vocabularies (fieldwork 1889, published 1897) are the earliest scholarly lists for the living languages, but which Arawan lects they cover was not verified (ARW-002, U). All other data are 20th–21st century (Chapman & Derbyshire 1991; Dixon 2004; Dienst 2014).
7. **Archaeology:** the upper-Purus geometric earthworks (Acre geoglyphs, c.2000–700 BP; Schaan et al. 2012, doi:10.1179/0093469012Z.00000000012) lie in or next to the Arawan/Purus-Arawakan/Panoan area. They have no linguistic attribution, and ARW-006 is listed only so it is not mistaken for a calibration.

### What this implies
* The independent indications put Proto-Arawan at **~1.2–3.5 ky**: ASJP ~1.7 ky; naive glottochronology 2.5–3.6 ky; Jolkesky's scenario ≤ ~3.2 ky; "shallow" from KAHD. These are not independent of lexical rate assumptions, so none can serve as a calibration for a lexical method.
* **Comparable families by ASJP similarity** (a proxy for depth under one method) are better anchored:

  | family | ASJP similarity | ASJP age (BP1950) |
  |---|---|---|
  | Tupí-Guaraní | 33.22 | 1500 |
  | Tacanan | 32.37 | 1540 |
  | Quechuan | 29.77 | 1667 |
  | Arawan | 28.88 | 1714 |
  | Panoan | 27.24 | 1803 |
  | Katukinan | 25.29 | 1915 |

  Of these, Tupí-Guaraní has a Bayesian posterior (1762 BP, 847–2740) and colonial tips, and Quechua has colonial tips and Inca-horizon bounds. So the most defensible route for (b) is **a cross-family rate prior estimated from calibrated shallow families with similar ASJP similarity and comparable sociolinguistic setting**:
  * Tupí-Guaraní (tips only; Ferraz Gerardi et al. 2023);
  * Arawakan subclades (Michael et al. 2023 calibrations);
  * Uto-Aztecan subclades (Greenhill et al. 2023);
  * Polynesian (short chronology);
  * Pama-Nyungan/Wati.

  Then report Arawan's root as a posterior and check it against the 1.2–3.5 ky window.
* **Neighbours:**
  * Katukinan: ASJP 1915 BP.
  * Harakmbut: no date found. Adelaar's Harakmbut–Katukinan link would be deep.
  * Panoan: ASJP 1803 BP. No Bayesian date. Zariquiey et al. 2023 does no dating.
  * Tacanan: ASJP 1540 BP.
  * Purus Arawakan (Apurinã, Yine, Iñapari): within Arawakan; Michael et al. 2023 has no Purus-specific calibration among the rows transcribed.
  * Tupian: ASJP 3535 BP; Tupí-Guaraní as above.
  * Chapacuran: ASJP 1881 BP.

  Arawakan → Arawan loans (Jolkesky 2016, Table 64) give only relative chronology. Dixon (2004:1) rejects an Arawakan–Arawan genetic relation.
* **Not found / still to check:**
  * any time-depth remark in the body of Dixon 2004 (only p.3 table seen via Dienst);
  * Kaufman 1990 per-family glottochronological depths (not seen; U);
  * Aikhenvald 1999 and Dixon 1999 (*The Amazonian Languages*) remarks;
  * any ASJP-tree (Müller et al.) node ages for Arawan subgroups;
  * a full KAHD release.

## 5. Overlaps and duplicates (deliberate, flagged)
* HOL rows repeat Holman's 15 in-scope calibration points with Holman's own rationale (which fork B could not read). They overlap NA-002/NA-021 (Inuit), NA-008 (Wakashan), NA-015 (MVS), MA-002 (Cholan), AN-022/023 (Chamic, Malayo-Chamic), AN-028/029 (Malagasy), AN-013–015 (East Polynesian), AN-003 (MP), AU-02 (Pama-Nyungan), SA-026 (TG coastal). **Use one row per node.**
* ASJP-037 (Austronesian) vs AN-001/036; ASJP-033 (Uto-Aztecan) vs NA-020; ASJP-034 (Mayan) vs MA-001: these pairs are method comparisons, useful as benchmark disagreement.

## 6. Main conflicts between studies
| node | study 1 | study 2 |
|---|---|---|
| Proto-East Polynesian | Gray 2009: 1150–1800 | Wilmshurst 2011 implies ~700–930 |
| Hawaii settlement | Wilmshurst: AD 1219–1266 | Athens 2014: AD 940–1130 |
| Proto-Oceanic | 3200–3600 (Gray 2009/2010) | 3000–3300 (Gray 2011) |
| PMP | 3600–4500 (Gray 2009) | N. Luzon pottery 5430–4290 cal BP |
| Chamic | Sidwell, Holman: 1550 (5th c. CE) | Gray 1800–2500; Maurits 2020 rejects it |
| Proto-Mayan | Kaufman c.4150 BP | ASJP 2170 BP |
| Thule/Inuit | ~1000 CE (classical) | ≥1200 CE (Friesen & Arnold 2008) |
| Tupí-Guaraní crown | Bayesian 1762 (847–2740) | Paraná pottery 2090–1740 cal BP; Rodrigues c.2500 |
| Pama-Nyungan root | Bayesian 5671 (Wati-dependent) | Evans & Jones / Holman 4450; ASJP 4245 |
| Arawan root | ASJP ~1714 | naive glottochronology 2.5–3.6 ky (derived) |
| Melanesian Pidgin | stabilised 1840s (Keesing 1988) | 1863+ labour trade (Mühlhäusler) |

## 7. Principal sources (DOIs verified via Crossref or the publisher where given)
* Holman, E.W. et al. 2011. Automated dating of the world's language families based on lexical similarity. *Current Anthropology* 52(6):841–875. doi:10.1086/662127 (full text read: calibration list, Tables 1, 3–7).
* Gray, R.D., Drummond, A.J. & Greenhill, S.J. 2009. *Science* 323:479–483. doi:10.1126/science.1166858.
* Greenhill, S.J., Drummond, A.J. & Gray, R.D. 2010. *PLoS ONE* 5:e9573. doi:10.1371/journal.pone.0009573.
* Greenhill, S.J. et al. 2017. *PNAS*. doi:10.1073/pnas.1700388114.
* Maurits, L. et al. 2020. *J. Language Evolution* 5:17–38. doi:10.1093/jole/lzz009.
* Wilmshurst, J.M. et al. 2011. *PNAS* 108:1815–1820. doi:10.1073/pnas.1015876108.
* Athens, J.S., Rieth, T.M. & Dye, T.S. 2014. *American Antiquity* 79:144–155. doi:10.7183/0002-7316.79.1.144.
* Burley, D., Weisler, M.I. & Zhao, J. 2012. *PLoS ONE* 7:e48769. doi:10.1371/journal.pone.0048769.
* Denham, T., Bronk Ramsey, C. & Specht, J. 2012. *Archaeology in Oceania* 47:39–46. doi:10.1002/j.1834-4453.2012.tb00113.x.
* Carson, M.T. 2020. *Radiocarbon* 62:1733–1754. doi:10.1017/RDC.2020.89.
* Cochrane, E.E. et al. 2021. *PLoS ONE* 16:e0251407. doi:10.1371/journal.pone.0251407.
* Crowther, A. et al. 2016. *PNAS* 113:6635–6640. doi:10.1073/pnas.1522714113.
* Barber, I.G. et al. 2016. *J. Archaeological Science: Reports* 10:96–109. doi:10.1016/j.jasrep.2016.08.040.
* Greenhill, S.J. et al. 2023. A recent northern origin for the Uto-Aztecan family. *Language* 99(1):81–107. doi:10.1353/lan.2023.0006.
* Friesen, T.M. & Arnold, C.D. 2008. *American Antiquity* 73(3). doi:10.1017/s0002731600046850.
* Houston, S., Robertson, J. & Stuart, D. 2000. *Current Anthropology* 41. doi:10.1086/300142.
* Kaufman, T. & Justeson, J. 2009. *Ancient Mesoamerica*. doi:10.1017/s0956536109990113.
* Campbell, L. & Kaufman, T. 1976. *American Antiquity* 41. doi:10.2307/279044.
* Justeson, J. & Kaufman, T. 1993. *Science* 259:1703. doi:10.1126/science.259.5102.1703.
* Michael, L. et al. 2023. Deriving calibrations for Arawakan using archaeological evidence. *Interface Focus* 13:20220049. doi:10.1098/rsfs.2022.0049.
* Walker, R.S. & Ribeiro, L.A. 2011. *Proc. R. Soc. B* 278:2562–2567. doi:10.1098/rspb.2010.2579 (undated).
* Ferraz Gerardi, F. et al. 2023. Lexical phylogenetics of the Tupí-Guaraní family. *PLOS ONE* 18(6):e0272226. doi:10.1371/journal.pone.0272226.
* Michael, L. et al. 2015. *LIAMES* 15(2):193–221. doi:10.20396/liames.v15i2.8642301 (undated).
* Bonomo, M. et al. 2015. *Quaternary International* 356:54–73. doi:10.1016/j.quaint.2014.10.050.
* Marsh, E.J. et al. 2017. *Radiocarbon* 59:117–140. doi:10.1017/RDC.2016.118.
* Heggarty, P. & Beresford-Jones, D. 2010. *Boletín de Arqueología PUCP* 14:29–60. doi:10.18800/boletindearqueologiapucp.201001.002.
* Rodrigues, A.D. 1964. *Revista de Antropologia* 12:99–104. doi:10.11606/2179-0892.ra.1964.110739.
* Bouckaert, R., Bowern, C. & Atkinson, Q.D. 2018. *Nature Ecology & Evolution* 2:741–749. doi:10.1038/s41559-018-0489-3.
* Bowern, C. 2012. *Proc. R. Soc. B*. doi:10.1098/rspb.2012.1842.
* Denham, T. et al. 2003. *Science* 301:189. doi:10.1126/science.1085255.
* Arawan:
  * Dixon, R.M.W. 2004. Proto-Arawá phonology. *Anthropological Linguistics* 46:1–83 (no DOI).
  * Dienst, S. 2008. *LIAMES* 8:61–67. doi:10.20396/liames.v8i1.1471.
  * Chandless, W. 1869. *J. Royal Geographical Society* 39:296–311. doi:10.2307/1798554.
  * Ferraz Gerardi, F., Aragon, C.C. & Reichert, S. 2022. *JOHD* 8:18. doi:10.5334/johd.80.
  * Jolkesky, M.P.V. 2016. PhD thesis, UnB (repositorio.unb.br/handle/10482/21671).
  * Schaan, D. et al. 2012. *J. Field Archaeology* 37:132–142. doi:10.1179/0093469012Z.00000000012.
  * Ehrenreich, P. 1897. *Materialien zur Sprachenkunde Brasiliens: Vokabulare von Purus-Stämmen*. Berlin: Asher.
* Colonial primary sources (tips) are cited in full in each TSV row.

## 8. Known weaknesses of this compilation
* The Gray 2009 SOM and the published Greenhill 2023 supplement were not read directly (secondary copy and preprint respectively).
* Grade-U rows (13) include:
  * San Vitores 1668 Chamorro;
  * Proto-Mixe-Zoque;
  * Proto-Aymaran;
  * Tukanoan;
  * Southern Jê;
  * Cumanagoto/Chaima;
  * Omagua-Kokama;
  * Kaurna 1840;
  * Ndyuka;
  * Negerhollands/Cape Verde dates;
  * ARW-002/004/006.
* Some glottocodes are left `.` where they could not be verified. Michael et al. 2023 codes are sometimes dialect-level, e.g. `hist1234` Historical Baure.
* Sampled-ancestor constraints (Michael et al. 2023) have no native clepsydra type. The TSV splits each into a crown soft maximum and a stem soft minimum, which loses information.
