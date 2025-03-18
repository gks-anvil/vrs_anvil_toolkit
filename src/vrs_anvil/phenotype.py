"""Utilize phenotype mapping and ontology traversal to generate dependency graphs
to support more flexible phenotype searches during cohort construction.
"""

from collections import defaultdict
from enum import Enum
from wags_tails import HpoData, MondoData
import fastobo


class OntologyPrefix(str, Enum):
    MONDO = "MONDO"
    HP = "HP"


class OntologyIndex:
    """Container for HPO and MONDO term lookups.

    >>> from vrs_anvil.phenotype import OntologyIndex
    >>> index = OntologyIndex()
    >>> index.get_child_terms("HP:0001618")
    {"HP:0001618", "HP:0007024", "HP:0012049"}

    """

    def __init__(self):
        """Construct ontology index.

        Ingest MONDO and HPO, build mappings between parent and child terms.
        """
        self.dependency_map = defaultdict(list)
        self.label_map = {}

        mondo_reader = self._get_mondo_reader()
        hpo_reader = self._get_hpo_reader()

        self._ingest_ontology(mondo_reader, OntologyPrefix.MONDO)
        self._ingest_ontology(hpo_reader, OntologyPrefix.HP)

    def _ingest_ontology(self, reader, prefix: OntologyPrefix):
        for term in reader:
            term_id = str(term.id)
            if not term_id.startswith(prefix.value):
                continue
            for clause in term:
                clause_tag = clause.raw_tag()
                if clause_tag == "is_a":
                    self.dependency_map[clause.raw_value()].append(term_id)
                elif clause_tag == "name":
                    self.label_map[clause.raw_value().upper()] = term_id

    @staticmethod
    def _get_mondo_reader():
        mondo_file = MondoData().get_latest()[0]
        return fastobo.iter(str(mondo_file.absolute()))

    @staticmethod
    def _get_hpo_reader():
        hpo_file = HpoData().get_latest()[0]
        return fastobo.iter(str(hpo_file.absolute()))

    def _get_dependency_set(self, parent_term) -> set[str]:
        children = {parent_term}
        for child in self.dependency_map[parent_term]:
            children |= self._get_dependency_set(child)
        return children


    def get_child_terms(self, term: str) -> set[str]:
        """Retrieve child terms for an ontology term, i.e. all terms which have an
        `is_a` relationship with the given term.

        Args:
            term: HPO or MONDO ontology term, including prefix

        Returns:
            set of child terms
        """
        if term not in self.dependency_map:
            raise KeyError

        return self._get_dependency_set(term)



class HpoOrganSystemCodes(str, Enum):
    GENITOURINARY = "HP:0000119"
    HEAD_NECK = "HP:0000152"
    EYE = "HP:0000478"
    EAR = "HP:0000598"
    NERVOUS = "HP:0000707"
    BREAST = "HP:0000769"
    ENDOCRINE = "HP:0000818"
    PRENATAL_BIRTH = "HP:0001197"
    GROWTH = "HP:0001507"
    INTEGUMENT = "HP:0001574"
    VOICE = "HP:0001608"
    CARDIOVASCULAR = "HP:0001626"
    BLOOD = "HP:0001871"
    METABOLISM_HOMEOSTASIS = "HP:0001939"
    RESPIRATORY = "HP:0002086"
    NEOPLASM = "HP:0002664"
    IMMUNE = "HP:0002715"
    DIGESTIVE = "HP:0025031"
    CONSTITUTIONAL = "HP:0025142"
    CELLULAR = "HP:0025354"
    MUSCULOSKELETAL = "HP:0033127"
    LIMBS = "HP:0040064"
    THORACIC = "HP:0045027"
