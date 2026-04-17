def circular_reaction_fingerprint(
    graph: CondensedReactionGraph,
    radius: int = 3,
    n_bits: int = 1024,
    count: bool = True,
    only_active_atoms: bool = False,
):
    gen = _reaction_circular_generator(graph)

    all_colors: list[np.ndarray] = []
    for colors, r in zip(gen, range(radius + 1)):
        if only_active_atoms:
            active_atoms = graph.active_atoms(additional_layer=r)

            colors = colors[np.fromiter(active_atoms, dtype=int)]
        all_colors.append(colors)

    fp = modulo_fold(np.concatenate(all_colors), n_bits=n_bits, count=count)
    return fp


def circular_stereo_reaction_fingerprint(
    graph: StereoCondensedReactionGraph,
    radius: int = 3,
    n_bits: int = 1024,
    count: bool = True,
    only_active_atoms: bool = False,
):
    gen = _stereo_reaction_circular_generator(graph)

    all_colors: list[np.ndarray] = []
    for colors, r in zip(gen, range(radius + 1)):
        if only_active_atoms:
            active_atoms = graph.active_atoms(additional_layer=r)

            colors = colors[np.fromiter(active_atoms, dtype=int)]
        all_colors.append(colors)

    fp = modulo_fold(np.concatenate(all_colors), n_bits=n_bits, count=count)
    return fp
