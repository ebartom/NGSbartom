from pkg_resources import resource_stream
import ConfigParser
#import eddlib

def load_parameters(default_param_filename='default_parameters.conf',
                    non_default_config_file=None):
    if non_default_config_file is not None:
        foo_config = non_default_config_file
    else:
        foo_config = resource_stream(__name__, default_param_filename)
    c = ConfigParser.ConfigParser()
    c.readfp(foo_config)
    d = {}
    d['fraq_ibins'] = c.getfloat('EDD config', 'required_fraction_of_informative_bins')
    d['ci_method'] = c.get('EDD config', 'p_hat_ci_method')
    d['ci_lim'] = c.getfloat('EDD config', 'max_CI_value')
    d['log_ratio_bin_size'] = c.getint('EDD config',
                                            'log_ratio_bin_size') * 1000
    assert 0.0 <= d['fraq_ibins'] <= 1.0
    assert 0.0 <= d['ci_lim'] <= 1.0
    return d
