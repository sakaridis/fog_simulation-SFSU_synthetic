function t = transmission_postprocessing_guided_filter(t_in, I, varargin)
%TRANSMISSION_POSTPROCESSING_GUIDED_FILTER  Process input transmission map with
%guided filter. Parameters for the filter are set to their recommended values.

window_size = 41;
mu = 1e-3;
t = clip_to_unit_range(transmission_guided_filtering(t_in, I, window_size, mu));

end

