from typing import List

import numpy as np
import sklearn.cluster
import distance
from loguru import logger
from spellchecker import SpellChecker

# CACHE & MANUALLY CURATED SYNONYM RESOLUTION
spellcheck_cache = {
    "h. sapiens": "homo sapiens",
    "howo sapiens": "homo sapiens",
    "human": "homo sapiens",
    "hunan": "homo sapiens",
    "homo sapiens": "homo sapiens",  # prevents spellchecker from correcting "homo sapiens" into "homosapien"
    "homosapien": "homo sapiens",
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


#   ############################### CORRECT USA REGIONS  #################################

USA_state_postal_codes = {
        'AL': 'Alabama',
        'AK': 'Alaska',
        'AZ': 'Arizona',
        'AR': 'Arkansas',
        'CA': 'California',
        'CO': 'Colorado',
        'CT': 'Connecticut',
        'DE': 'Delaware',
        'FL': 'Florida',
        'GA': 'Georgia',
        'HI': 'Hawaii',
        'ID': 'Idaho',
        'IL': 'Illinois',
        'IN': 'Indiana',
        'IA': 'Iowa',
        'KS': 'Kansas',
        'KY': 'Kentucky',
        'LA': 'Louisiana',
        'ME': 'Maine',
        'MD': 'Maryland',
        'MA': 'Massachusetts',
        'MI': 'Michigan',
        'MN': 'Minnesota',
        'MS': 'Mississippi',
        'MO': 'Missouri',
        'MT': 'Montana',
        'NE': 'Nebraska',
        'NV': 'Nevada',
        'NH': 'New Hampshire',
        'NJ': 'New Jersey',
        'NM': 'New Mexico',
        'NY': 'New York',
        'NC': 'North Carolina',
        'ND': 'North Dakota',
        'OH': 'Ohio',
        'OK': 'Oklahoma',
        'OR': 'Oregon',
        'PA': 'Pennsylvania',
        'RI': 'Rhode Island',
        'SC': 'South Carolina',
        'SD': 'South Dakota',
        'TN': 'Tennessee',
        'TX': 'Texas',
        'UT': 'Utah',
        'VT': 'Vermont',
        'VA': 'Virginia',
        'WA': 'Washington',
        'WV': 'West Virginia',
        'WI': 'Wisconsin',
        'WY': 'Wyoming',
        'AS': 'American Samoa',
        'DC': 'District of Columbia',
        'FM': 'Federated States of Micronesia',
        'GU': 'Guam',
        'MH': 'Marshall Islands',
        'MP': 'Northern Mariana Islands',
        'PW': 'Palau',
        'PR': 'Puerto Rico',
        'VI': 'Virgin Islands'
}
USA_state_names_upper_case = {v.upper(): v for v in USA_state_postal_codes.values()}    # still a dictionary


def correct_usa_regions(region: str):
    if not region or 'unknown' in region.lower():
        return None
    region_parts_upper_case = [x.strip().upper() for x in region.split(',', maxsplit=1)]   # REGION PARTS ARE UPPER CASE

    if len(region_parts_upper_case) == 1:  # no comma in region
        # it may be a postal code
        state_name = USA_state_postal_codes.get(region_parts_upper_case[0])

        # it could be one of the special cases
        if state_name is None:
            state_name = {
                # special cases
                'CALIFORNI': 'California',
                'SLIDELL LA': 'Lousiana',
            }.get(region_parts_upper_case[0])
        return state_name or region

    elif len(region_parts_upper_case) == 2:   # one comma
        # check left hand side of comma

        # one of the parts could be the state name
        state_name = USA_state_names_upper_case.get(region_parts_upper_case[0])
        if not state_name:
            state_name = USA_state_names_upper_case.get(region_parts_upper_case[1])

        # one of the parts could be the postal code of the state
        if not state_name:
            state_name = USA_state_postal_codes.get(region_parts_upper_case[0])
        if not state_name:
            state_name = USA_state_postal_codes.get(region_parts_upper_case[1])
        if not state_name:
            state_name = region_parts_upper_case[0].capitalize()
            if '/' in state_name:
                try:
                    state_name = state_name[:state_name.rindex('/') - 1].rstrip()
                except:
                    pass
        return state_name

    else:
        logger.warning(f"Correction of USA country names. Region '{region}' is not handled")
        return region

