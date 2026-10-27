"""The annotations of the recurring compartments of the models.

An annotation is a list of `(qualifier, resource)` pairs: the qualifier is a
biological qualifier of `sbmlutils.metadata.BQB` and the resource a term of an
ontology, written as `collection/term` of identifiers.org. sbmlutils writes
them as the RDF of the element they are given to, e.g.
`Compartment("Vli", ..., annotations=compartments["li"])` for the liver.

The PBPK library (`pkpdsim.pbpk.annotation`) takes the terms of its organs,
blood pools, lumen and tubule segments and sampling sites from here and adds
the terms of the PBPK ontology.
"""

from sbmlutils.metadata import BQB

compartments: dict[str, list[tuple[BQB, str]]] = {
    # the extracellular fluid
    "fluid volume": [
        (BQB.IS, "mesh/D045604"),  # extracellular fluid
        (BQB.IS_VERSION_OF, "ncit/C13236"),  # body fluid or substance
    ],
    # liver
    "li": [
        (BQB.IS, "fma/FMA:7197"),  # liver
        (BQB.IS, "bto/BTO:0000759"),  # liver
        (BQB.IS, "ncit/C12392"),  # liver
        (BQB.IS, "uberon/UBERON:0002107"),  # liver
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
    # rest of the body: the tissues of the organism which no other
    # compartment of the model holds, a part of the organism
    "re": [
        (BQB.IS_VERSION_OF, "fma/FMA:62955"),  # anatomical entity
        (BQB.IS_PART_OF, "ncit/C14250"),  # Organism
        (BQB.IS_PART_OF, "uberon/UBERON:0000468"),  # multicellular organism
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
        (BQB.IS, "uberon/UBERON:0000948"),  # heart
    ],
    # kidneys
    "ki": [
        (BQB.IS, "fma/FMA:7203"),  # kidney
        (BQB.IS, "bto/BTO:0000671"),  # kidney
        (BQB.IS, "ncit/C12415"),  # kidney
        (BQB.IS, "uberon/UBERON:0002113"),  # kidney
    ],
    # lung
    "lu": [
        (BQB.IS, "fma/FMA:7195"),  # lung
        (BQB.IS, "bto/BTO:0000763"),  # lung
        (BQB.IS, "ncit/C12468"),  # lung
        (BQB.IS, "uberon/UBERON:0002048"),  # lung
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
    # rest vein: the veins of the rest of the body
    "rev": [
        (BQB.IS, "bto/BTO:0000131"),  # plasma
        (BQB.IS, "ncit/C13356"),  # plasma
        (BQB.IS_PART_OF, "uberon/UBERON:0004582"),  # venous system
        (BQB.IS_PART_OF, "ncit/C12814"),  # vein
    ],
    # kidney vein
    "kiv": [
        (BQB.IS, "bto/BTO:0000131"),  # plasma
        (BQB.IS, "ncit/C13356"),  # plasma
        (BQB.IS_PART_OF, "uberon/UBERON:0001140"),  # renal vein
        (BQB.IS_PART_OF, "ncit/C33462"),  # renal vein
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
        (BQB.IS, "uberon/UBERON:0001134"),  # skeletal muscle tissue
    ],
    # forearm
    "fo": [
        (BQB.IS, "ncit/C32628"),  # forearm
        (BQB.IS, "fma/FMA:9663"),  # forearm
    ],
    # forearm vein: the superficial veins of the forearm, the sampling site
    "fov": [
        (BQB.IS, "bto/BTO:0000131"),  # plasma
        (BQB.IS, "ncit/C13356"),  # plasma
        (BQB.IS_PART_OF, "ncit/C12814"),  # vein
        (BQB.IS_PART_OF, "ncit/C32628"),  # forearm
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
        (BQB.IS, "uberon/UBERON:0001969"),  # blood plasma
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
    # the organs of the PBPK models which have no key above
    "br": [
        (BQB.IS, "uberon/UBERON:0000955"),  # brain
        (BQB.IS, "ncit/C12439"),  # brain
    ],
    "si": [
        (BQB.IS, "uberon/UBERON:0002108"),  # small intestine
        (BQB.IS, "ncit/C12386"),  # small intestine
    ],
    "la": [
        (BQB.IS, "uberon/UBERON:0000059"),  # large intestine
        (BQB.IS, "ncit/C12379"),  # large intestine
    ],
    "sk": [
        (BQB.IS, "uberon/UBERON:0002097"),  # skin of body
        (BQB.IS, "ncit/C12470"),  # skin
    ],
    "ad": [
        (BQB.IS, "uberon/UBERON:0001013"),  # adipose tissue
        (BQB.IS, "ncit/C12472"),  # adipose tissue
    ],
    "bo": [
        (BQB.IS, "uberon/UBERON:0002481"),  # bone tissue
        (BQB.IS, "ncit/C12366"),  # bone
    ],
    "go": [
        (BQB.IS, "uberon/UBERON:0000991"),  # gonad
        (BQB.IS, "ncit/C12725"),  # gonad
    ],
    # the segments of the large intestine
    "caecum": [
        (BQB.IS, "uberon/UBERON:0001153"),  # caecum
        (BQB.IS, "ncit/C12381"),  # cecum
    ],
    "ascending colon": [
        (BQB.IS, "uberon/UBERON:0001156"),  # ascending colon
        (BQB.IS, "ncit/C12265"),  # ascending colon
    ],
    "transverse colon": [
        (BQB.IS, "uberon/UBERON:0001157"),  # transverse colon
        (BQB.IS, "ncit/C12385"),  # transverse colon
    ],
    "descending colon": [
        (BQB.IS, "uberon/UBERON:0001158"),  # descending colon
        (BQB.IS, "ncit/C12268"),  # descending colon
    ],
    "sigmoid colon": [
        (BQB.IS, "uberon/UBERON:0001159"),  # sigmoid colon
        (BQB.IS, "ncit/C12384"),  # sigmoid colon
    ],
    "rectum": [
        (BQB.IS, "uberon/UBERON:0001052"),  # rectum
        (BQB.IS, "ncit/C12390"),  # rectum
    ],
    # the segments of the renal tubule
    "proximal tubule": [
        (BQB.IS, "uberon/UBERON:0004134"),  # proximal tubule
    ],
    "loop of Henle": [
        (BQB.IS, "uberon/UBERON:0001288"),  # loop of Henle
    ],
    "distal tubule": [
        (BQB.IS, "uberon/UBERON:0004135"),  # distal tubule
    ],
    "collecting duct": [
        (BQB.IS, "uberon/UBERON:0001232"),  # collecting duct of renal tubule
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
