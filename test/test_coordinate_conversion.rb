require 'rubygems'
require 'test/unit'
require 'lib/geokit'

class TestCoordinateConversion < Test::Unit::TestCase
  def setup
    @delta = 0.00005
  end

  def test_to_cartesian
    r = Geokit::LatLng.send(:units_sphere_multiplier, {})
    point = Geokit::LatLng.new(80, 30)
    coords = point.to_cartesian

    expected = [0.8528685 * r, 0.4924039 * r, 0.1736482 * r]
    expected.each_with_index do |p, i|
      assert_in_delta p, coords[i], @delta
    end
  end

  def test_to_polar
    r = Geokit::LatLng.send(:units_sphere_multiplier, {})
    point = [0.8528685 * r, 0.4924039 * r, 0.1736482 * r]
    coords = Geokit::LatLng.send(:to_polar, *point)

    assert_in_delta 80, coords.first, @delta
    assert_in_delta 30, coords.last, @delta
  end
end
