python scripts/setup_db.py
python scripts/plot_quality_control.py
python scripts/propagate_flags.py
python scripts/homogenise_giraffe.py
python scripts/plot_giraffe_homogenisation_diagnostics.py
python scripts/homogenise_uves.py
python scripts/plot_uves_homogenisation_diagnostics.py
python scripts/plot_homogenisation_quality_control.py
tar -cvzf outputs/ges-corot-homogenisation-figures.tar.gz figures/*
python scripts/shipit.py
