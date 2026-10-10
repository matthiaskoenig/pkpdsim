"""The annotations of the recurring compartments of the models.

An annotation is a list of `(qualifier, resource)` pairs: the qualifier is a
biological qualifier of `sbmlutils.metadata.BQB` and the resource a term of an
ontology, written as `collection/term` of identifiers.org. sbmlutils writes
them as the RDF of the element they are given to, e.g.
`Compartment("Vli", ..., annotations=compartments["li"])` for the liver.
"""

from sbmlutils.metadata import BQB

compartments: dict[str, list[tuple[BQB, str]]] = {
    "fluid volume": [
        # FIXME: annotations
    ],
    # liver
    "li": [
        (BQB.IS, "fma/FMA:7197"),  # liver
        (BQB.IS, "bto/BTO:0000759"),  # liver
        (BQB.IS, "ncit/C12392"),  # liver
    ],
    "lyso": [
        (BQB.IS, "fma/FMA:63836"),  # lysosome
        (BQB.IS, "ncit/C13253"),  # lysosome
    ],
    # basolateral membrane
    "basolateral": [
        (BQB.IS, "fma/FMA:84669"),  # Basolateral plasma membrane
        (BQB.IS, "go/GO:0016323"),  # basolateral plasma membrane
    ],
    # apical membrane
    "apical": [
        (BQB.IS, "fma/FMA:84666"),  # Apical plasma membrane
        (BQB.IS, "go/GO:0016324"),  # apical plasma membrane
    ],
    # bile
    "bi": [
        (BQB.IS, "fma/FMA:62971"),  # portion of bile
        (BQB.IS, "ncit/C13192"),  # bile
    ],
    # rest
    # FIXME: rest compartment is not the organism
    "re": [
        (BQB.IS_VERSION_OF, "fma/FMA:62955"),  # anatomical entity
        (BQB.IS_PART_OF, "ncit/C14250"),  # Organism
    ],
    # gastrointestinal tract
    "gi": [
        (BQB.IS, "fma/FMA:71132"),  # gastrointestinal tract
        (BQB.IS, "bto/BTO:0000511"),  # gastrointestinal tract
        (BQB.IS, "ncit/C34082"),  # gastrointestinal tract
    ],
    # gut/intestine
    "gu": [
        (BQB.IS, "bto/BTO:0000545"),  # gut
        (BQB.IS, "ncit/C12736"),  # intestine
        (BQB.IS, "fma/FMA:7199"),  # intestine
        (BQB.IS, "bto/BTO:0000648"),  # intestine
    ],
    "gu_lumen": [
        (BQB.IS, "fma/FMA:14586"),  # Lumen of intestine
        (BQB.IS, "uberon/UBERON:0018543"),  # lumen of intestine
    ],
    "duodenum": [
        (BQB.IS, "fma/FMA:7206"),  # duodenum
        (BQB.IS, "uberon/UBERON:0002114"),  # duodenum
    ],
    "jejunum": [
        (BQB.IS, "fma/FMA:7207"),  # jejunum
        (BQB.IS, "uberon/UBERON:0002115"),  # jejunum
    ],
    "ileum": [
        (BQB.IS, "fma/FMA:7208"),  # ileum
        (BQB.IS, "uberon/UBERON:0002116"),  # ileum
    ],
    "colon": [
        (BQB.IS, "fma/FMA:14543"),  # colon
        (BQB.IS, "uberon/UBERON:0001155"),  # colon
    ],
    # heart
    "he": [
        (BQB.IS, "ncit/C12727"),  # heart
        (BQB.IS, "fma/FMA:7088"),  # heart
        (BQB.IS, "bto/BTO:0000562"),  # heart
    ],
    # kidneys
    "ki": [
        (BQB.IS, "fma/FMA:7203"),  # kidney
        (BQB.IS, "bto/BTO:0000671"),  # kidney
        (BQB.IS, "ncit/C12415"),  # kidney
    ],
    # lung
    "lu": [
        (BQB.IS, "fma/FMA:7195"),  # lung
        (BQB.IS, "bto/BTO:0000763"),  # lung
        (BQB.IS, "ncit/C12468"),  # lung
    ],
    "ve": [
        (BQB.IS, "bto/BTO:0000131"),  # plasma
        (BQB.IS, "ncit/C13356"),  # plasma
        (BQB.IS_PART_OF, "fma/FMA:50723"),  # venous tree
    ],
    "ar": [
        (BQB.IS, "bto/BTO:0000131"),  # plasma
        (BQB.IS, "ncit/C13356"),  # plasma
        (BQB.IS_PART_OF, "fma/FMA:50720"),  # arterial tree
    ],
    # portal vein
    "po": [
        (BQB.IS, "bto/BTO:0000131"),  # plasma
        (BQB.IS, "ncit/C13356"),  # plasma
        (BQB.IS_PART_OF, "fma/FMA:66645"),  # portal vein
    ],
    # hepatic vein
    "hv": [
        (BQB.IS, "bto/BTO:0000131"),  # plasma
        (BQB.IS, "ncit/C13356"),  # plasma
        (BQB.IS_PART_OF, "fma/FMA:14337"),  # hepatic vein
    ],
    # renal vein
    "rv": [
        (BQB.IS, "bto/BTO:0000131"),  # plasma
        (BQB.IS, "ncit/C13356"),  # plasma
        (BQB.IS_PART_OF, "bto/BTO:0002681"),  # renal vein
        (BQB.IS_PART_OF, "ncit/C33462"),  # renal vein
        (BQB.IS_PART_OF, "fma/FMA:14334"),  # renal vein
    ],
    # rest vein
    "rev": [
        (BQB.IS, "bto/BTO:0000131"),  # plasma
        (BQB.IS, "ncit/C13356"),  # plasma
        # fixme: be more specific
    ],
    # kidney vein
    "kiv": [
        (BQB.IS, "bto/BTO:0000131"),  # plasma
        (BQB.IS, "ncit/C13356"),  # plasma
        # fixme: be more specific
    ],
    "stomach": [
        (BQB.IS, "fma/FMA:7148"),  # stomach
        (BQB.IS, "ncit/C12391"),  # stomach
        (BQB.IS, "uberon/UBERON:0000945"),  # stomach
        (BQB.IS, "bto/BTO:0001307"),  # stomach
    ],
    # spleen
    "sp": [
        (BQB.IS, "fma/FMA:7196"),  # spleen
        (BQB.IS, "ncit/C12432"),  # spleen
        (BQB.IS, "uberon/UBERON:0002106"),  # spleen
        (BQB.IS, "bto/BTO:0001281"),  # spleen
    ],
    # muscle
    "mu": [
        (BQB.IS, "fma/FMA:30316"),  # muscle
        (BQB.IS, "bto/BTO:0000887"),  # muscle
        (BQB.IS, "ncit/C13056"),  # muscle
    ],
    # forearm
    "fo": [
        (BQB.IS, "ncit/C32628"),  # forearm
        (BQB.IS, "fma/FMA:9663"),  # forearm
    ],
    # forearm vein
    "fov": [
        (BQB.IS, "bto/BTO:0000131"),  # plasma
        (BQB.IS, "ncit/C13356"),  # plasma
        # fixme: be more specific
    ],
    # pancreas
    "pa": [
        (BQB.IS, "fma/FMA:7198"),  # pancreas
        (BQB.IS, "ncit/C12393"),  # pancreas
        (BQB.IS, "uberon/UBERON:0001264"),  # pancreas
        (BQB.IS, "bto/BTO:0000988"),  # pancreas
    ],
    "feces": [
        (BQB.IS, "fma/FMA:64183"),  # feces
        (BQB.IS, "ncit/C13234"),  # feces
        (BQB.IS, "bto/BTO:0000440"),  # feces
    ],
    "plasma": [
        (BQB.IS, "ncit/C13356"),  # plasma
        (BQB.IS, "bto/BTO:0000131"),  # plasma
    ],
    "blood": [
        (BQB.IS, "fma/FMA:9670"),  # portion of blood
        (BQB.IS, "ncit/C12434"),  # blood
        (BQB.IS, "bto/BTO:0000089"),  # blood
    ],
    "urine": [
        (BQB.IS, "fma/FMA:12274"),  # portion of urine
        (BQB.IS, "ncit/C13283"),  # urine
        (BQB.IS, "bto/BTO:0001419"),  # urine
    ],
    "parenchyma": [
        (BQB.IS, "fma/FMA:45732"),  # parenchyma
        (BQB.IS, "ncit/C74601"),  # parenchyma
        (BQB.IS, "bto/BTO:0001539"),  # parenchyma
    ],
    "parietal cell": [
        (BQB.IS, "ncit/C12594"),  # parietal cell
        (BQB.IS, "snomedct/57041003"),  # gastric parietal cell
    ],
    "gastric acid": [
        (BQB.IS, "fma/FMA:62972"),  # portion of gastric juice
        (BQB.IS, "ncit/C94192"),  # gastric acid
        (BQB.IS, "omit/0006944"),  # gastric acid
    ],
    "plasma membrane": [
        (BQB.IS, "ncit/C13735"),  # Plasma membrane
        (BQB.IS, "fma/FMA:63841"),  # Plasma membrane
        (BQB.IS, "go/GO:0005886"),  # plasma membrane
    ],
}
"""The annotations of the compartments, keyed by their abbreviation in the models.

`li` liver, `ki` kidney, `gu` gut, `ve` and `ar` the venous and the arterial
blood, `po` the portal vein, `urine`, `feces`, ...
"""
