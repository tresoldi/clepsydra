"""
Concept ID mapping for the reference builder.

Maps external concept identifiers (Concepticon IDs, Swadesh list labels)
to Clepsydra canonical concept IDs. The C inference engine uses only the
canonical ID; external mappings live in this separate table.

NOTE: In a real build, this table is populated from Concepticon. For the
synthetic fixture, we generate a small Swadesh-style concept list.
"""

# Canonical concept IDs are stable ASCII identifiers prefixed "C".
# A small subset for synthetic testing; the real build uses the full
# Swadesh-200 / Concepticon mapping.

# (canonical_id, label, concepticon_id, word_class, status)
SYNTHETIC_CONCEPTS = [
    ("C000948", "WATER", "948", "noun", "active"),
    ("C000149", "TWO", "149", "numeral", "active"),
    ("C000025", "FIRE", "25", "noun", "active"),
    ("C000031", "HAND", "31", "noun", "active"),
    ("C000049", "DOG", "49", "noun", "active"),
    ("C000011", "ONE", "11", "numeral", "active"),
    ("C000029", "EYE", "29", "noun", "active"),
    ("C000058", "FISH", "58", "noun", "active"),
    ("C000069", "SUN", "69", "noun", "active"),
    ("C000071", "TONGUE", "71", "noun", "active"),
    ("C000075", "TOOTH", "75", "noun", "active"),
    ("C000079", "BLOOD", "79", "noun", "active"),
    ("C000093", "BONE", "93", "noun", "active"),
    ("C000096", "HORN", "96", "noun", "active"),
    ("C000100", "BREAST", "100", "noun", "active"),
    ("C000104", "HEART", "104", "noun", "active"),
    ("C000118", "MOUNTAIN", "118", "noun", "active"),
    ("C000125", "STAR", "125", "noun", "active"),
    ("C000130", "STONE", "130", "noun", "active"),
    ("C000140", "WIND", "140", "noun", "active"),
    ("C000141", "DIRTY", "141", "adjective", "active"),
    ("C000145", "SMALL", "145", "adjective", "active"),
    ("C000149b", "THREE", "150", "numeral", "active"),
    ("C000151", "PERSON", "151", "noun", "active"),
    ("C000159", "TREE", "159", "noun", "active"),
    ("C000165", "BIG", "165", "adjective", "active"),
    ("C000169", "COME", "169", "verb", "active"),
    ("C000175", "DIE", "175", "verb", "active"),
    ("C000180", "EAT", "180", "verb", "active"),
    ("C000185", "GO", "185", "verb", "active"),
    ("C000190", "HEAR", "190", "verb", "active"),
    ("C000195", "I", "195", "pronoun", "active"),
    ("C000200", "YOU", "200", "pronoun", "active"),
    ("C000205", "WHO", "205", "pronoun", "active"),
    ("C000210", "WHAT", "210", "pronoun", "active"),
    ("C000215", "NOT", "215", "adverb", "active"),
    ("C000220", "ALL", "220", "adverb", "active"),
    ("C000225", "AND", "225", "conjunction", "active"),
    ("C000230", "IN", "230", "preposition", "active"),
    ("C000235", "WITH", "235", "preposition", "active"),
    ("C_CUSTOM_001", "CUSTOM_CONCEPT", ".", "unknown", "active"),
]


def get_concepts():
    """Return the concept list as a list of dicts."""
    return [
        {
            "concept_id": cid,
            "label": label,
            "external_namespace": "Concepticon" if cnc != "." else ".",
            "external_id": cnc,
            "word_class": wc,
            "status": status,
        }
        for (cid, label, cnc, wc, status) in SYNTHETIC_CONCEPTS
    ]


def canonical_ids():
    return [c[0] for c in SYNTHETIC_CONCEPTS]