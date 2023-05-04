mod paths;

pub fn get_all_paths(
    contours: &pxu::Contours,
    consts: pxu::CouplingConstants,
) -> Vec<pxu::path::SavedPath> {
    paths::ALL_PATHS
        .iter()
        .map(|f| f(&contours, consts))
        .collect::<Vec<_>>()
}
