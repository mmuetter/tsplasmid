def resolve_significance_table(M, parent_group = [], level = 0):
    sig_groups = []
    for i, row in M.iterrows():
        row = row.dropna()
        parent_group, group = resolve_row(row, M, level, parent_group = parent_group[:level])
        sig_groups.extend(group)
    return sig_groups


def resolve_row(row, M, level, parent_group = []):
  #  print("\nRow:", row.name, "\n Parent: ", parent_group)
    row = row.dropna()
    groups = row.index
    M_sub = M.copy()
    not_significant =  row.values == False
    M_sub = M_sub.loc[row.index, row.index]
    M_sub = M_sub.loc[not_significant, not_significant]

    if len(M_sub) == 0:
   #     print("return:",set(parent_group + [row.name]))
        return parent_group, [set(parent_group + [row.name])]
    elif len(M_sub) == 1:
    #    print("return:",set(parent_group + [row.name] + list(M_sub.index)))
        return parent_group, [set(parent_group + [row.name] + list(M_sub.index)) ]
    else:
        parent_group.append(row.name)
        return parent_group, resolve_significance_table(M_sub, parent_group = parent_group, level =  level+1)
    
