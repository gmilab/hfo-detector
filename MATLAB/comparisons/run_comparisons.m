% Donos data
hfo_dset_donos = hfodat.load_donos();

donos_rl_confusion = hfo_dset_donos.run_all_ripplelab();
donos_rl_summary = hfodat.summarize_confusion(donos_rl_confusion);

donos_zrh_confusion = hfo_dset_donos.compute_confusion(hfo_dset_donos.run_zurich());
donos_zrh_summary = hfodat.summarize_confusion(donos_zrh_confusion);

donos_hsc_confusion = hfo_dset_donos.compute_confusion(hfo_dset_donos.run_hscdetector(), 't_threshold', 0.5);
donos_hsc_summary = hfodat.summarize_confusion(donos_hsc_confusion);

% Anywave data
hfo_dset_anywave = hfodat.load_anywave();

anywave_rl_confusion = hfo_dset_anywave.run_all_ripplelab();
anywave_rl_summary = hfodat.summarize_confusion(anywave_rl_confusion);

anywave_zrh_confusion = hfo_dset_anywave.compute_confusion(hfo_dset_anywave.run_zurich());
anywave_zrh_summary = hfodat.summarize_confusion(anywave_zrh_confusion);

anywave_hsc_confusion = hfo_dset_anywave.compute_confusion(hfo_dset_anywave.run_hscdetector(), 't_threshold', 0.5);
anywave_hsc_summary = hfodat.summarize_confusion(anywave_hsc_confusion);


% format percentages
donos_rl_summary.PPV_fmt = arrayfun(@(x) sprintf('%.1f', x * 100), donos_rl_summary.PPV, 'UniformOutput', false);
donos_rl_summary.SEN_fmt = arrayfun(@(x) sprintf('%.1f', x * 100), donos_rl_summary.SEN, 'UniformOutput', false);
donos_zrh_summary.PPV_fmt = arrayfun(@(x) sprintf('%.1f', x * 100), donos_zrh_summary.PPV, 'UniformOutput', false);
donos_zrh_summary.SEN_fmt = arrayfun(@(x) sprintf('%.1f', x * 100), donos_zrh_summary.SEN, 'UniformOutput', false);
donos_hsc_summary.PPV_fmt = arrayfun(@(x) sprintf('%.1f', x * 100), donos_hsc_summary.PPV, 'UniformOutput', false);
donos_hsc_summary.SEN_fmt = arrayfun(@(x) sprintf('%.1f', x * 100), donos_hsc_summary.SEN, 'UniformOutput', false);

anywave_rl_summary.PPV_fmt = arrayfun(@(x) sprintf('%.1f', x * 100), anywave_rl_summary.PPV, 'UniformOutput', false);
anywave_rl_summary.SEN_fmt = arrayfun(@(x) sprintf('%.1f', x * 100), anywave_rl_summary.SEN, 'UniformOutput', false);
anywave_zrh_summary.PPV_fmt = arrayfun(@(x) sprintf('%.1f', x * 100), anywave_zrh_summary.PPV, 'UniformOutput', false);
anywave_zrh_summary.SEN_fmt = arrayfun(@(x) sprintf('%.1f', x * 100), anywave_zrh_summary.SEN, 'UniformOutput', false);
anywave_hsc_summary.PPV_fmt = arrayfun(@(x) sprintf('%.1f', x * 100), anywave_hsc_summary.PPV, 'UniformOutput', false);
anywave_hsc_summary.SEN_fmt = arrayfun(@(x) sprintf('%.1f', x * 100), anywave_hsc_summary.SEN, 'UniformOutput', false);

% Display output
display(donos_rl_summary)
display(donos_zrh_summary)
display(donos_hsc_summary)
display(anywave_rl_summary)
display(anywave_zrh_summary)
display(anywave_hsc_summary)


