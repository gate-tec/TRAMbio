from typing import List
from TRAMbio.util.structure_library.pdb_entry_map import EntryMap


def parse_entry_by_dict(entry: str, mapping: List[EntryMap]):
    result = {}
    for i, ele in enumerate(mapping):
        if ele['line'][0] < len(entry) and (len(ele['line']) == 1 or ele['line'][1] <= len(entry)):
            if len(ele['line']) == 1:
                value = entry[ele['line'][0]:].strip()
            else:
                value = entry[ele['line'][0]: ele['line'][1]].strip()
            try:
                result[ele['id']] = ele['type'](value)
            except ValueError:
                # parsing error
                if ele['required']:
                    return None
                else:
                    result[ele['id']] = None
        elif ele['required']:
            return None
        else:
            result[ele['id']] = None

    return result
