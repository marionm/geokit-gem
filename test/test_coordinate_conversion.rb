require 'rubygems'
require 'test/unit'
require 'lib/geokit'

class TestCoordinateConversions < Test::Unit::TestCase
  def test_to_cartesian
    radius = Geokit::Mappable::EARTH_RADIUS_IN_MILES
    { [0, 0]   => [radius, 0, 0],
      [0, 180] => [-radius, 0, 0],
      [90, 0]  => [0, 0, radius],
      [-90, 0] => [0, 0, -radius],
      [0, 90]  => [0, radius, 0],
      [0, -90] => [0, -radius, 0],
      [45, 0]  => [radius * Math.sin(Math::PI * 0.25), 0, radius * Math.sin(Math::PI * 0.25)]
    }.each do |input, expected|
      actual = Geokit::LatLng.new(*input).to_cartesian
      expected.each_with_index do |coord, i|
        assert_in_delta coord, actual[i], 0.00000005
      end
    end
  end

  def test_to_polar
    radius = Geokit::Mappable::EARTH_RADIUS_IN_MILES
    { [radius, 0, 0]  => [0, 0],
      [-radius, 0, 0] => [0, 180],
      [0, radius, 0]  => [0, 90],
      [0, -radius, 0] => [0, -90],
      [0, 0, radius]  => [90, 0],
      [0, 0, -radius] => [-90, 0],
      [radius * Math.sin(Math::PI * 0.25), 0, radius * Math.sin(Math::PI * 0.25)] => [45, 0],
      [-radius * Math.sin(Math::PI * 0.25), radius * Math.sin(Math::PI * 0.25), 0] => [0, 135]
    }.each do |input, expected|
      actual = Geokit::LatLng.send(:to_polar, *input)
      assert_in_delta expected.first, actual.first, 0.00000005
      assert_in_delta expected.last,  actual.last,  0.00000005
    end
  end
end
