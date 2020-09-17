from typing import List

import numpy as np
import sklearn.cluster
import distance
from spellchecker import SpellChecker

# CACHE & MANUALLY CURATED SYNONYM RESOLUTION
spellcheck_cache = {
    "howo sapiens": "homo sapiens",
    "human": "homo sapiens",
    "ae. aegypti": "aedes aegypti",
    "aedes aegypti mosquitoes": "aedes aegypti",
    "aedes albopictus, female": "aedes albopictus",
    "culicidae: culex sp.": "culex sp.",
    "culicidae: culex vaxus": "culex vaxus",
    "enviromental": "environment",
    "pipistrellus cf. hesperidus": "pipistrellus hesperidus",
    "mosquitoes": "mosquito"
}
spell = SpellChecker()


#   ################   CLUSTERING
def cluster(input_words):
    words = np.asarray(input_words) #So that indexing with a list will work
    lev_similarity = -1*np.array([[distance.levenshtein(w1,w2) for w1 in words] for w2 in words])

    affprop = sklearn.cluster.AffinityPropagation(affinity="precomputed", damping=0.5, random_state=0, convergence_iter=100)
    affprop.fit(lev_similarity)
    clusters = []
    for cluster_id in np.unique(affprop.labels_):
        exemplar = words[affprop.cluster_centers_indices_[cluster_id]]
        cluster = np.unique(words[np.nonzero(affprop.labels_==cluster_id)])
        clusters.append((exemplar, cluster))
    return clusters


def print_clusters(clusters):
    for c in clusters:
        exemplar = c[0]
        cluster = c[1]
        cluster_str = ", ".join(cluster)
        print(" - *%s:* %s" % (exemplar, cluster_str))


#   #############      TYPOS CORRECTION
def correct_typos_in_list(input_words: List[str]) -> List[str]:
    global spellcheck_cache

    # find those words that may be misspelled
    misspelled = spell.unknown(input_words)

    result = []
    for word in input_words:
        # Get the one `most likely` answer
        correction = spellcheck_cache.get(word.lower())
        if not correction:
            correction = spell.correction(word).lower()
            spellcheck_cache[word] = correction
        result.append(correction)
    return result


def correct_typos(input_word) -> str:
    global spellcheck_cache
    input_word = input_word.lower()
    correction = spellcheck_cache.get(input_word)
    if not correction:
        correction = spell.correction(input_word).lower()
        spellcheck_cache[input_word] = correction
    return correction


def _just_a_test():
    my_words = [
        "Macaca fascicularis",
        "Rhinolophus sp.",
        "field Aedes albopictus",
        "dromedary camel",
        "Hypsugo savii",
        "Homo sapien",
        "dromedary",
        "Ae. aegypti",
        "mosquitoes",
        "Culex sp.",
        "Aedes albopictus",
        "swine",
        "camelus dromedarius",
        "Chaerephon pumilus",
        "Aedes albopictus, female",
        "dengue fever patient",
        "Camel",
        "Macaca mulatta",
        "arthropod",
        "gorilla",
        "mosquito suspension",
        "Homo sapiense",
        "dengue",
        "Homo sapines",
        "goat",
        "Rhogeessa tumida",
        "Sus scrofa",
        "Vero cell",
        "bat",
        "canine",
        "chimpanzee",
        "domestic donkey",
        "Lama glama",
        "Pythium insidiosum",
        "monkey",
        "Camelus dromedarius",
        "mouse",
        "Environment",
        "Carollia perspicillata",
        "Mops condylurus",
        "Molossus pretiosus",
        "Homo sapiens",
        "cynomolgus macaque",
        "palm civet",
        "Culex quinquefasciatus",
        "mosquito",
        "Marmosa murina",
        "Culicidae: Culex vaxus",
        "pig",
        "howo sapiens",
        "Phyllostomus discolor",
        "camel",
        "Molossus sinaloae",
        "Proechimys cuvieri",
        "Aedes aegypti",
        "Bovidae",
        "Aedes aegypti mosquitoes",
        "Canis familiaris",
        "Mustela lutreola",
        "Mus musculus",
        "Rhinolophus ferrumequinum",
        "sentinel monkey",
        "Felis catus",
        "Pipistrellus kuhlii",
        "guinea pig",
        "Neoromicia capensis",
        "Pipistrellus cf. hesperidus",
        "Enviromental",
        "Culicidae: Culex sp.",
        "Panthera tigris jacksoni",
        "Molossus rufus",
        "Glossophaga soricina",
        "Homo sapience",
        "primate",
        "Human",
        "Didelphis marsupialis"
    ]
    print(f"words in input: {len(my_words)}")
    # correct_words = correct_typos(my_words)
    # print(f"words in correct_words: {len(correct_words)}")

    from datetime import datetime
    total_time_start = datetime.now()
    for f in my_words:
        print(f"{f} ->\t\t", end='')
        start = datetime.now()
        correct = correct_typos(f)
        print(f"{correct}\t\ttime: {datetime.now()-start}")
    total_time_end = datetime.now()
    print(f'total time: {total_time_end- total_time_start}\t\taverage time: {(total_time_end-total_time_start)/len(my_words)}')

