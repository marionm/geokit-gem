require 'test/unit'
require 'rubygems'
require 'lib/geokit'

# TODO: Better and more tests would be nice
class PolygonalMethodsTest < Test::Unit::TestCase
  def setup
    #Points generated by http://www.geocodezip.com/map-markers_ConvexHull_Polygon.asp
    @locations = [
      [37.461371,-122.154351],
      [37.458874,-122.131191],
      [37.455747,-122.126103],
      [37.422689,-122.107946],
      [37.435285,-122.132227],
      [37.423376,-122.139946],
      [37.422976,-122.146517],
      [37.439795,-122.151437],
      [37.441023,-122.181853],
      [37.452271,-122.180637]
    ].map { |coords| Geokit::LatLng.new(*coords) }
  end

  def test_contrived_polygon_example
    polygon = Geokit::LatLng.polygon(@locations.shuffle)
    assert_equal @locations, polygon
  end

  def test_contrived_convex_hull_example
    expected = [
      @locations[0], @locations[1], @locations[2],
      @locations[3], @locations[6], @locations[8],
      @locations[9]
    ]
    hull = Geokit::LatLng.convex_hull(@locations.shuffle)
    assert_equal expected, hull
  end
end

class TestPolygonClockwise < Test::Unit::TestCase
  def test_clockwise
    [
      [[0, 0], [0, 1], [-1, -2]], [[0, 0], [0, 1], [-1, 0]], [[0, 0], [0, 1], [-1, 1]], [[0, 0], [0, 1], [-1, 2]],
      [[0, 0], [0, -1], [1, -2]], [[0, 0], [0, -1], [1, -1]], [[0, 0], [0, -1], [1, 0]], [[0, 0], [0, -1], [1, 1]],
      [[30, 54], [67, 54], [-20, 55]], [[40, 10], [10, 10], [30, 0]], [[40, 50], [60, 20], [64, 30]]
    ].each do |points|
      assert Geokit::LatLng.send(:clockwise?, *points.map { |p| Geokit::LatLng.new(*p) })
    end
  end

  def test_counter_clockwise
    [
      [[0, 0], [0, 1], [1, -2]], [[0, 0], [0, 1], [1, 0]], [[0, 0], [0, 1], [1, 1]], [[0, 0], [0, 1], [1, 2]],
      [[0, 0], [0, -1], [-1, -2]], [[0, 0], [0, -1], [-1, -1]], [[0, 0], [0, -1], [-1, 0]], [[0, 0], [0, -1], [-1, 1]],
      [[30, 54], [67, 54], [-20, 52]], [[40, 10], [10, 10], [30, 85]],
      [[40, 50], [60, 20], [60, 10]]
    ].each do |points|
      assert !Geokit::LatLng.send(:clockwise?, *points.map { |p| Geokit::LatLng.new(*p) })
    end
  end
end

class TestCentroid < Test::Unit::TestCase
  def test_centroid_of_one_point
    point = Geokit::LatLng.new(65, 43)
    assert_equal point, Geokit::LatLng.centroid([point])
  end

  def test_centroid_of_two_points
    points = [
      a = Geokit::LatLng.new(65, 43),
      b = Geokit::LatLng.new(31, -18)
    ]
    midpoint = Geokit::LatLng.midpoint_between(a, b)
    actual = Geokit::LatLng.centroid(points)

    assert_in_delta midpoint.lat, actual.lat, 0.000005
    assert_in_delta midpoint.lng, actual.lng, 0.000005
  end

  # TODO: More tests would be nice
  def test_centroid_of_three_points
    points = [
      Geokit::LatLng.new(0, -45),
      Geokit::LatLng.new(0, 45),
      Geokit::LatLng.new(90, 0)
    ]
    expected = Geokit::LatLng.new(35.26438, 0)
    actual = Geokit::LatLng.centroid(points)
    assert_in_delta expected.lat, actual.lat, 0.00005
    assert_in_delta expected.lng, actual.lng, 0.00005

    points = [
      Geokit::LatLng.new(0, -45),
      Geokit::LatLng.new(0, 45),
      Geokit::LatLng.new(-90, 0)
    ]
    expected = Geokit::LatLng.new(-35.26438, 0)
    actual = Geokit::LatLng.centroid(points)
    assert_in_delta expected.lat, actual.lat, 0.00005
    assert_in_delta expected.lng, actual.lng, 0.00005
  end

  def test_centroid_of_x_points_as_polygon
    points = [
      Geokit::LatLng.new(0, -45),
      Geokit::LatLng.new(0, 45),
      Geokit::LatLng.new(89, 1),
      Geokit::LatLng.new(89, -1)
    ]
    expected = Geokit::LatLng.new(34.9, 0)
    actual = Geokit::LatLng.centroid(points, :method => :polygon)
    assert_in_delta expected.lat, actual.lat, 0.05
    assert_in_delta expected.lng, actual.lng, 0.0005

    points = [
      Geokit::LatLng.new(0, -2),
      Geokit::LatLng.new(-1, 0),
      Geokit::LatLng.new(0, 2),
      Geokit::LatLng.new(-2, 2),
      Geokit::LatLng.new(-3, 0),
      Geokit::LatLng.new(-2, -2)
    ]
    expected = Geokit::LatLng.new(-1.3335, 0)
    actual = Geokit::LatLng.centroid(points, :method => :polygon)
    assert_in_delta expected.lat, actual.lat, 0.005
    #FIXME: This test correctly fails - there is a bug in anchor selections for the polygon mehod
    assert_in_delta expected.lng, actual.lng, 0.005
  end

  def test_centroid_of_x_points_as_convex_hull
    points = [
      Geokit::LatLng.new(0, -2),
      Geokit::LatLng.new(-1, 0),
      Geokit::LatLng.new(0, 2),
      Geokit::LatLng.new(-2, 2),
      Geokit::LatLng.new(-3, 0),
      Geokit::LatLng.new(-2, -2)
    ]

    expected = Geokit::LatLng.new(-1.2668, 0)
    actual = Geokit::LatLng.centroid(points, :method => :convex_hull)

    assert_in_delta expected.lat, actual.lat, 0.005
    assert_in_delta expected.lng, actual.lng, 0.0005
  end
end

class TestArea < Test::Unit::TestCase
  def test_greater_but_close_to_planar_approx
    a = [0, 0]
    b = [10, 0]
    c = [0, 10]

    approx = ((10 * Geokit::Mappable::MILES_PER_LATITUDE_DEGREE) ** 2) / 2
    actual = Geokit::LatLng.area([a, b, c])

    assert approx < actual
    assert_in_delta approx, actual, 1800
  end

  def test_misleading_headings
    #The heading difference between b to a and b to c is ~296, which is not the angle at b
    a = [5, 0]
    b = [-3, 0]
    c = [-2, -2]

    lat_miles = Geokit::Mappable::MILES_PER_LATITUDE_DEGREE
    approx = 8 * lat_miles * 2 * lat_miles / 2
    actual = Geokit::LatLng.area([a, b, c])

    assert approx < actual
    assert_in_delta approx, actual, 101
  end
end